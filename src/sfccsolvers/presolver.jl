abstract type AbstractPreSolverCache end

"""
    SFCCPreSolver{TCache} <: SFCCSolver

A fast SFCC "solver" which pre-builds an interpolant for the freeze curve in terms of enthalpy, θ(H).
"""
struct SFCCPreSolver{TCache} <: SFCCSolver
    cache::TCache
    Tmin::Float64
    errtol::Float64
    reinit::Bool # if true, always reinitialize the presolver cache when initialize! is called
    SFCCPreSolver(cache, Tmin, errtol, reinit=true) = new{typeof(cache)}(cache, Tmin, errtol, reinit)
    """
        SFCCPreSolver(cache::AbstractPreSolverCache; Tmin=-60.0, errtol=1e-4)

    Constructs a new `SFCCPreSolver` with minimum temperature `Tmin` and tolerance `errtol`.
    Enthalpy values below `H(Tmin)` under the given freeze curve will be extrapolated with a
    constant/flat function. `errtol` determines the permitted local error in the interpolant.
    """
    function SFCCPreSolver(cache::AbstractPreSolverCache=SFCCPreSolverCacheND(); Tmin=-60.0, errtol=1e-4, reinit=true)
        new{typeof(cache)}(cache, Tmin, errtol, reinit)
    end
end

"""
    SFCCPreSolverCache1D{TF} <: AbstractPreSolverCache

Efficient "presolver" implementation that builds a 1D interpolant for `θw(H)` over an adaptive range of `H` values determined based on local approximation error.
Note that this presolver implementation is **only valid when all other freeze curve parameters are held constant** and will produce incorrect results otherwise.
"""
mutable struct SFCCPreSolverCache1D{TF} <: AbstractPreSolverCache
    f_H::TF # θw(H) interpolant
    initialized::Bool
    function SFCCPreSolverCache1D()
        # initialize with dummy functions to get type information;
        # this is just to make sure that the compiler can still infer all of the types.
        x  = -3e8:1e6:3e8
        dummy_f = _build_interpolant1D(x, zeros(length(x)))
        return new{typeof(dummy_f)}(dummy_f, false)
    end
end
function _build_interpolant1D(H, θw)
    return Interpolations.extrapolate(
        Interpolations.interpolate((Vector(H),), θw, Interpolations.Gridded(Interpolations.Linear())),
        Interpolations.Flat()
    )
end
"""
    initialize!(solver::SFCCPreSolver{<:SFCCPreSolverCache1D}, fc::SFCC, hc::F; fc_kwargs...)

Initializes this `SFCCPreSolver` by building interpolants for `fc` and its derivatives. `θtot` and `θsat` are assumed fixed.
`hc` must be heat capacity as a function of liquid water content.
"""
function initialize!(solver::SFCCPreSolver{<:SFCCPreSolverCache1D}, fc::SFCC, hc::F; sat=1.0, θsat=SoilWaterVolume(fc).θsat, fc_kwargs...) where {F}
    if solver.cache.initialized && !solver.reinit
        # skip re-initialization
        return nothing
    end
    # pre-solve freeze curve;
    # note that this is only valid given that the following assumptions hold:
    # 1) none of the freeze curve parameters (e.g. soil properties) change
    # 2) soil properties are uniform in `soil`
    let Tmin = solver.Tmin,
        Tmax = SoilFreezeThawProperties(fc).Tₘ,
        Lsl = SoilFreezeThawProperties(fc).Lsl,
        ρw = SoilWaterVolume(fc).ρw,
        L = ρw*Lsl,
        f_kwargs = (; θsat, fc_kwargs...),
        f(T) = fc(T; f_kwargs...),
        Hmin = FreezeCurves.enthalpy(Tmin, hc(f(Tmin), sat*θsat, θsat), L, f(Tmin)),
        Hmax = FreezeCurves.enthalpy(Tmax, hc(f(Tmax), sat*θsat, θsat), L, f(Tmax));
        # residual as a function of T and H
        resid(T,H) = FreezeCurves.temperature_residual(fc, hc, L, H, T, sat; f_kwargs...)
        function solve(H,T₀)
            local T = nlsolve(T -> resid(T[1],H)[1], [T₀]).zero[1]
            local θw = f(T)
            return (; T, θw)
        end
        function deriv(T) # implicit partial derivative w.r.t H as a function of T
            local θw, ∂θ∂T = ∇(f, T)
            # get ∂H∂T
            _, ∂H∂T = ∇(T -> let θw=f(T); FreezeCurves.enthalpy(T, hc(θw, sat*θsat, θsat), L, θw); end, T)
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
        ∂θw∂T₀, ∂θw∂H₀ = deriv(T[1])
        ∂θw∂H = [∂θw∂H₀]
        @assert isfinite(H[1]) && isfinite(θw[1]) "H=$H, θw=$θw"
        while H[end] < Hmax
            # find the optimal step size
            ϵ = Inf
            ΔH = L/10 # initially set to (relatively) large value, 10% of latent heat of fusion
            while abs(ϵ) > solver.errtol
                ϵ = step(ΔH, H[end], θw[end], ∂θw∂H[end], T[end])
                # iteratively halve the step size until error tolerance is satisfied
                ΔH *= 0.5
            end
            Hnew = H[end] + ΔH
            @assert isfinite(Hnew) "isfinite(ΔH) failed; H=$(H[end]), T=$(T[end]), ΔH=$ΔH"
            opt = solve(Hnew, T[end])
            ∂θw∂Hᵢ, ∂θw∂Tᵢ = deriv(opt.T)
            push!(H, Hnew)
            push!(θw, opt.θw)
            push!(T, opt.T)
            push!(∂θw∂H, ∂θw∂Hᵢ)
        end
        solver.cache.f_H = _build_interpolant1D(H, θw)
        solver.cache.initialized = true
    end
    return nothing
end

# SFCC solve implementation for presolver
function sfccsolve(
    obj::SFCCInverseEnthalpyObjective,
    solver::SFCCPreSolver{<:SFCCPreSolverCache1D},
    ::Any,
    ::Val{return_all}=Val{true}()
) where {return_all}
    @assert solver.cache.initialized "solver not yet initialized"
    θsat = get(obj.f_kwargs, :θsat, SoilWaterVolume(obj.f))
    H = @SVector[obj.H]
    H_dual = dual(H, typeof(obj))[1]
    θw_dual = solver.cache.f_H(H_dual)
    C_dual = obj.hc(θw_dual, obj.sat*θsat, θsat)
    T_dual = (H_dual - obj.L*θw_dual) / C_dual
    ∂θw∂H = ForwardDiff.partials(θw_dual)[1]
    ∂T∂H = ForwardDiff.partials(T_dual)[1]
    ∂θw∂T = ∂θw∂H / (∂T∂H + eps())
    T = ForwardDiff.value(T_dual)
    θw = ForwardDiff.value(θw_dual)
    C = ForwardDiff.value(C_dual)
    return return_all ? (; T, θw, C, ∂θw∂T, ∂T∂H, ∂θw∂H) : T
end

"""
    SFCCPreSolverCacheND{N,TF} <: AbstractPreSolverCache

More general "presolver" implementation that builds a ND lookup-table/interpolant for `θw(H, sat, args...)` where `args` is
a set of all other freeze curve parameters passed to `initialize!`. Note that this implementation currently does not perform
any error control so there approximation error is not guaranteed.
"""
mutable struct SFCCPreSolverCacheND{N,TF} <: AbstractPreSolverCache
    lut::TF # T(H) interpolant
    initialized::Bool
    function SFCCPreSolverCacheND(; H=-1e8:1e5:3.34e8, sat=0.00:0.01:1.0, θsat=0.1:0.1:1.0, other_coords...)
        dummy_f(; kwargs...) = NaN
        lut = build_lut(dummy_f, :H => H, :sat => sat, :θsat => θsat, other_coords...)
        N = length(lut.coords)
        return new{N,typeof(lut)}(lut, false)
    end
end

function initialize!(solver::SFCCPreSolver{<:SFCCPreSolverCacheND{N}}, fc::SFCC, hc::F; sat=1.0, θsat=0.5, fc_kwargs...) where {N,F}
    if solver.cache.initialized && !solver.reinit
        # skip re-initialization
        return nothing
    end
    let Tₘ = SoilFreezeThawProperties(fc).Tₘ,
        Lsl = SoilFreezeThawProperties(fc).Lsl,
        ρw = SoilWaterVolume(fc).ρw,
        L = ρw*Lsl,
        f(T, sat; kwargs...) = fc(T, sat; kwargs...);
        # residual as a function of T and H
        resid(T, H, sat, f_kwargs) = FreezeCurves.temperature_residual(fc, hc, L, H, T, sat; f_kwargs...)
        function solve_T(;H, sat, θsat, kwargs...)
            f_kwargs = (; θsat, kwargs...)
            obj = SFCCInverseEnthalpyObjective(fc, f_kwargs, hc, L, H, sat)
            res = sfccsolve(obj, SFCCNewtonSolver(), Tₘ-1.0)
            return res.T
        end
        lut = build_lut(
            solve_T,
            :H => solver.cache.lut.coords.H,
            :sat => solver.cache.lut.coords.sat,
            :θsat => solver.cache.lut.coords.θsat,
            fc_kwargs...
        )
        solver.cache.lut = lut
        solver.cache.initialized = true
        return nothing
    end
end

function sfccsolve(
    obj::SFCCInverseEnthalpyObjective,
    solver::SFCCPreSolver{<:SFCCPreSolverCacheND},
    ::Any,
    ::Val{return_all}=Val{true}()
) where {return_all}
    @assert solver.cache.initialized "solver not yet initialized"
    θsat = get(obj.f_kwargs, :θsat, SoilWaterVolume(obj.f).θsat)
    θtot = obj.sat*θsat
    L = obj.L
    H = @SVector[obj.H]
    H_dual = dual(H, typeof(obj))[1]
    T_dual = IfElse.ifelse(
        obj.H > L*θtot,
        (H_dual - L*θtot) / obj.hc(θtot, θtot, θsat),
        solver.cache.lut(; H=H_dual, sat=obj.sat, obj.f_kwargs...)
    )
    θw_dual = obj.f(T_dual, obj.sat; obj.f_kwargs...)
    C_dual = obj.hc(θw_dual, obj.sat*θsat, θsat)
    ∂θw∂H = ForwardDiff.partials(θw_dual)[1]
    ∂T∂H = ForwardDiff.partials(T_dual)[1]
    ∂θw∂T = ∂θw∂H / (∂T∂H + eps())
    T = ForwardDiff.value(T_dual)
    θw = ForwardDiff.value(θw_dual)
    C = ForwardDiff.value(C_dual)
    return return_all ? (; T, θw, C, ∂θw∂T, ∂T∂H, ∂θw∂H) : T
end
