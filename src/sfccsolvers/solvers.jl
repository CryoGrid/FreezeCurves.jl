using ForwardDiff
using IfElse
using Interpolations
using NLsolve
using StaticArrays
using Requires
using Unitful

export sfccsolve, SFCCInverseEnthalpyObjective

"""
    AbstractSFCCObjective

Base type for non-linear SFCC optimization objectives.
"""
abstract type AbstractSFCCObjective end

"""
    SFCCInverseEnthalpyObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH,Tsat} <: AbstractSFCCObjective

Optimization objective for finding a temperature, `T`, wich resolves the conservation equation:
```math
0 = (H - Lθ)/C - T
```
given a fixed enthalpy value `H` and freeze curve arguments `f_args`.
"""
Base.@kwdef struct SFCCInverseEnthalpyObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH,Tsat} <: AbstractSFCCObjective
    f::TF
    f_kwargs::Tkwargs
    hc::Thc
    L::TL
    H::TH
    sat::Tsat = 1.0
end
@inline (obj::SFCCInverseEnthalpyObjective)(T) = first(FreezeCurves.temperature_residual(obj.f, obj.hc, obj.L, adstrip(obj.H), T, adstrip(obj.sat); obj.f_kwargs...))

"""
    initialize!(solver::SFCCSolver, fc::SFCC, hc; fc_kwargs...)

Initializes `solver` (if necessary) for the freeze curve function `fc` with arguments `fc_kwargs` and heat capacity function `hc`.
Default implementation does nothing.
"""
initialize!(solver::SFCCSolver, fc::SFCC, hc; fc_kwargs...) = nothing

"""
    sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀, ::Val{return_all}=Val{true}()) where {return_all}

Solve the given objective `obj` using `solver` and initial guess `x₀`. If `return_all=true`, then `sfccsolve` should return a named tuple
with at least the temperature solution `T`, the liquid water content `θw`, the heat capacity `C`, and the liquid water content deriviative
`∂θw∂T` defined. Solver-specific additional values may also be included.
"""
sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀, ::Val{return_all}=Val{true}()) where {return_all} = error("not implemented for objective $(typeof(obj)) with solver $(typeof(solver))")

"""
adstrip(x::Number)
adstrip(x::ForwardDiff.Dual)

adstrip extracts the underlying numeric value from `x` if `x` is a `ForwardDiff.Dual` number.
If `x` is a non-tracked numeric type, then `adstrip` simply returns `x`.
"""
adstrip(x::Number) = x
adstrip(x::ForwardDiff.Dual) = adstrip(ForwardDiff.value(x))
adstrip(::Nothing) = nothing

include("heatcap.jl")
export SFCCNewtonSolver
include("newton.jl")
export LUT, build_lut
include("lut.jl")
export SFCCPreSolver
include("presolver.jl")
