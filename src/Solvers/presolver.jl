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
    function SFCCPreSolver(; Tmin=-60.0, errtol=1e-4)
        cache = SFCCPreSolverCache1D()
        new{typeof(cache)}(cache, Tmin, errtol)
    end
end
mutable struct SFCCPreSolverCache1D{TFH,T∇FH,T∇FT}
    f_H::TFH # θw(H) interpolant
    dθwdH::T∇FH # derivative of θw(H) interpolant
    dθwdT::T∇FT # ∂θ∂T interpolant
    initialized::Bool
    function SFCCPreSolverCache1D()
        # initialize with dummy functions to get type information;
        # this is just to make sure that the compiler can still infer all of the types.
        x  = -3e8:1e6:3e8
        dummy_f = _build_interpolant(x, zeros(length(x)))
        dummy_∇f = first ∘ _interpgrad(dummy_f)
        return new{typeof(dummy_f),typeof(dummy_∇f),typeof(dummy_f)}(dummy_f, dummy_∇f, dummy_f, false)
    end
end
_interpgrad(f) = (args...) -> Interpolations.gradient(f, args...)
function _build_interpolant(H, θw)
    return Interpolations.extrapolate(
        Interpolations.interpolate((Vector(H),), θw, Interpolations.Gridded(Interpolations.Linear())),
        Interpolations.Flat()
    )
end
"""
    initialize!(solver::SFCCPreSolver{<:SFCCPreSolverCache1D}, fc::SFCCFunction, hc::F; fc_kwargs...)

Initializes this `SFCCPreSolver` by building interpolants for `fc` and its derivatives. `θtot` and `θsat` are assumed fixed.
`hc` must be heat capacity as a function of liquid water content.
"""
function initialize!(solver::SFCCPreSolver{<:SFCCPreSolverCache1D}, fc::SFCCFunction, hc::F; fc_kwargs...) where {F}
    # pre-solve freeze curve;
    # note that this is only valid given that the following assumptions hold:
    # 1) none of the freeze curve parameters (e.g. soil properties) change
    # 2) soil properties are uniform in `soil`
    let Tmin = solver.Tmin,
        Tmax = 0.0,
        ρw = SoilWaterProperties(fc).ρw,
        Lsl = SoilFreezeThawProperties(fc).Lsl,
        L = ρw*Lsl,
        f_kwargs = (; fc_kwargs...),
        f(T) = fc(T; f_kwargs...),
        Hmin = FreezeCurves.enthalpy(Tmin, hc(f(Tmin)), L, f(Tmin)),
        Hmax = FreezeCurves.enthalpy(Tmax, hc(f(Tmax)), L, f(Tmax));
        # residual as a function of T and H
        resid(T,H) = FreezeCurves.temperature_residual(fc, f_kwargs, hc, L, H, T)
        function solve(H,T₀)
            local T = nlsolve(T -> resid(T[1],H)[1], [T₀]).zero[1]
            local θw = f(T)
            return (; T, θw)
        end
        function deriv(T) # implicit partial derivative w.r.t H as a function of T
            local θw, ∂θ∂T = ∇(f, T)
            # get dHdT
            _, ∂H∂T = ∇(T -> let θw=f(T); FreezeCurves.enthalpy(T, hc(θw), L, θw); end, T)
            # valid by chain rule and inverse function theorem
            return ∂θ∂T / ∂H∂T, ∂θ∂T
        end
        function step(ΔH, H, θw, ∂θ∂H, T₀)
            # get first order estimate
            θest = θw + ΔH*∂θ∂H
            # get true θ at H + ΔH
            θsol = solve(H + ΔH, T₀).θw
            err = abs(θsol - θest)
            # return residual of error with target error
            return err
        end
        T = [Tmin]
        H = [Hmin]
        θw = [f(T[1])]
        dθwdT₀, dθwdH₀ = deriv(T[1])
        dθwdT = [dθwdT₀]
        dθwdH = [dθwdH₀]
        @assert isfinite(H[1]) && isfinite(θw[1]) "H=$H, θw=$θw"
        while H[end] < Hmax
            # find the optimal step size
            ϵ = Inf
            ΔH = L/10 # initially set to (relatively) large value, 10% of latent heat of fusion
            while abs(ϵ) > solver.errtol
                ϵ = step(ΔH, H[end], θw[end], dθwdH[end], T[end])
                # iteratively halve the step size until error tolerance is satisfied
                ΔH *= 0.5
            end
            Hnew = H[end] + ΔH
            @assert isfinite(Hnew) "isfinite(ΔH) failed; H=$(H[end]), T=$(T[end]), ΔH=$ΔH"
            opt = solve(Hnew, T[end])
            dθwdHᵢ, dθwdTᵢ = deriv(opt.T)
            push!(H, Hnew)
            push!(θw, opt.θw)
            push!(T, opt.T)
            push!(dθwdT, dθwdTᵢ)
            push!(dθwdH, dθwdHᵢ)
        end
        solver.cache.f_H = _build_interpolant(H, θw)
        solver.cache.dθwdH = first ∘ _interpgrad(solver.cache.f_H)
        solver.cache.dθwdT = _build_interpolant(T, dθwdT)
        solver.cache.initialized = true
    end
end
function sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCPreSolver, ::Any, ::Val{return_all}=Val{true}()) where {return_all}
    @assert solver.cache.initialized "solver not yet initialized"
    θw = solver.cache.f_H(obj.H)
    C = obj.hc(θw)
    T = (obj.H - obj.L*θw) / C
    dθwdH = solver.cache.dθwdH(obj.H)
    dθwdT = solver.cache.dθwdT(T)
    return return_all ? (; T, θw, C, dθwdT, dθwdH) : T
end
