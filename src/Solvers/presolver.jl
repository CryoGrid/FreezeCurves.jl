"""
    SFCCPreSolver{TCache} <: SFCCSolver

A fast SFCC "solver" which pre-builds an interpolant for the freeze curve in terms of enthalpy, θ(H).
Note that this solver is **only valid when all freeze curve parameters are held constant** and will
produce incorrect results otherwise.
"""
struct SFCCPreSolver{TCache} <: SFCCSolver
    cache::TCache
    Tmin::Float64
    errtol::Float64
    SFCCPreSolver(cache, Tmin, errtol) = new{typeof(cache)}(cache, Tmin, errtol)
    """
        SFCCPreSolver(Tmin=-60.0, errtol=1e-4)

    Constructs a new `SFCCPreSolver` with minimum temperature `Tmin` and integration step `dH`.
    Enthalpy values below `H(Tmin)` under the given freeze curve will be extrapolated with a
    constant/flat function. `errtol` determines the permitted local error in the interpolant.
    """
    function SFCCPreSolver(;Tmin=-60.0, errtol=1e-4)
        cache = SFCCPreSolverCache()
        new{typeof(cache)}(cache, Tmin, errtol)
    end
end
mutable struct SFCCPreSolverCache{TFH,T∇FH,T∇FT}
    f_H::TFH # θw(H) interpolant
    ∇f_H::T∇FH # derivative of θw(H) interpolant
    ∇f_T::T∇FT # ∂θ∂T interpolant
    function SFCCPreSolverCache()
        # initialize with dummy functions to get type information;
        # this is just to make sure that the compiler can still infer all of the types.
        x  = -3e8:1e6:3e8
        dummy_f = _build_interpolant(x, zeros(length(x)))
        dummy_∇f = first ∘ _interpgrad(dummy_f)
        return new{typeof(dummy_f),typeof(dummy_∇f),typeof(dummy_f)}(dummy_f, dummy_∇f, dummy_f)
    end
end
_interpgrad(f) = (args...) -> Interpolations.gradient(f, args...)
function _build_interpolant(Hs, θs)
    return Interpolations.extrapolate(
        Interpolations.interpolate((Vector(Hs),), θs, Interpolations.Gridded(Interpolations.Linear())),
        Interpolations.Flat()
    )
end
function initialize!(obj::SFCCTemperatureObjective, solver::SFCCPreSolver)
    # pre-solve freeze curve;
    # note that this is only valid given that the following assumptions hold:
    # 1) none of the freeze curve parameters (e.g. soil properties) change
    # 2) soil properties are uniform in `soil`
    let Tmin = solver.Tmin,
        Tmax = 0.0,
        f(T) = fc(T, θwi, θsat, f_args...),
        Hmin = enthalpy(Tmin, heatcap(f(Tmin)), L, f(Tmin)),
        Hmax = enthalpy(Tmax, heatcap(f(Tmax)), L, f(Tmax));
        # residual as a function of T and H
        resid(T,H) = temperature_residual(f, args, C, L, H, T)
        function solve(H,T₀)
            local T = nlsolve(T -> resid(T[1],H)[1], [T₀]).zero[1]
            θw = f(T)
            return (; T, θw)
        end
        function deriv(T) # implicit partial derivative w.r.t H as a function of T
            θw, ∂θ∂T = ∇(f, T)
            # get dHdT
            _, ∂H∂T = ∇(T -> let θw=f(T); enthalpy(T, heatcap(θw), L, θw); end, T)
            # valid by chain rule and inverse function theorem
            return ∂θ∂T / ∂H∂T, ∂θ∂T
        end
        function step(ΔH, H, θ, ∂θ∂H, T₀)
            # get first order estimate
            θest = θ + ΔH*∂θ∂H
            # get true θ at H + ΔH
            θsol = solve(H + ΔH, T₀).θw
            err = abs(θsol - θest)
            # return residual of error with target error
            return err
        end
        T = [Tmin]
        H = [Hmin]
        θ = [f(T[1])]
        dθdT₀, dθdH₀ = deriv(T[1])
        dθdT = [dθdT₀]
        dθdH = [dθdH₀]
        @assert isfinite(H[1]) && isfinite(θ[1]) "H=$H, θ=$θ"
        while H[end] < Hmax
            # find the optimal step size
            ϵ = Inf
            ΔH = L*θtot/10 # initially set to large value
            while abs(ϵ) > sfcc.solver.errtol
                ϵ = step(ΔH, H[end], θ[end], dθdH[end], T[end])
                # iteratively halve the step size until error tolerance is satisfied
                ΔH *= 0.5
            end
            Hnew = H[end] + ΔH
            @assert isfinite(Hnew) "isfinite(ΔH) failed; H=$(H[end]), T=$(T[end]), ΔH=$ΔH"
            opt = solve(Hnew, T[end])
            dθdHᵢ, dθdTᵢ = deriv(opt.T)
            push!(H, Hnew)
            push!(θ, opt.θw)
            push!(T, opt.T)
            push!(dθdT, dθdTᵢ)
            push!(dθdH, dθdHᵢ)
        end
        sfcc.solver.cache.f_H = _build_interpolant(H, θ)
        sfcc.solver.cache.∇f_H = first ∘ _interpgrad(sfcc.solver.cache.f_H)
        sfcc.solver.cache.∇f_T = _build_interpolant(T, dθdT)
    end
end
