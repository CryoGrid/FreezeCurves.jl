module FreezeCurves

using ForwardDiff
using IfElse
using Interpolations
using IntervalSets
using RecipesBase
using Reexport
using Requires
using StaticArrays
@reexport using Setfield: @set, @set!
@reexport using Unitful

import Flatten
import Unitful: ustrip

export FreezeCurve

function __init__()
    @require Turing="fce5fe82-541a-59a6-adf8-730c64b5f9a0" begin
        using .Turing, .Turing.Distributions
        export SFCCModel, sfccpriors
        include("inference/inference.jl")
    end
    # require NonlinearSolve.jl for generic nonlinear solver
    @require NonlinearSolve="8913a72c-1f9b-4ce2-8d82-65094dcecaec" begin
        using .NonlinearSolve
        export SFCCNonlinearSolver
        include("sfccsolvers/nonlinearsolve.jl")
    end
end

# convenience constants for temperature unit/quantity types
const TemperatureUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TemperatureQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TemperatureUnit}

"""
    FreezeCurve

Base type for freezing characteristic curves which relate temperature or enthalpy to unfrozen water content.
"""
abstract type FreezeCurve end

"""
    enthalpy(T, C, L, θw)

Computes the volumetric enthalpy [J/m^3] given T, C, L, θw.
"""
enthalpy(T, C, L, θw) = T*C + L*θw
"""
    normalize_temperature(x)
    normalize_temperature(x::TemperatureQuantity)

Converts temperature `x` to Kelvin. If `x` has units, `uconvert` is used. Otherwise, if `x` a general numeric type, it is assumed that `x` is in celsius.
"""
normalize_temperature(x) = x + 273.15
normalize_temperature(x::TemperatureQuantity) = uconvert(u"K", x)

"""
    heaviside(x)

Differentiable implementation of heaviside step function, i.e:

``h(x) = \\begin{cases} 1 & x ≥ 0 \\\\ 0 & x < 0 \\end{cases}``
"""
heaviside(x) = IfElse.ifelse(x >= zero(x), 1.0, 0.0)

include("math.jl")

export FreeWater
include("freewater.jl")

export SWRC, SoilWaterVolume, BrooksCorey, VanGenuchten
export inflectionpoint
include("swrc.jl")

export SFCC, SFCCSolver, SoilFreezeThawProperties
export PainterKarra, DallAmico, DallAmicoSalt, McKenzie, Westermann, Langer, Hu2020, PowerLaw
export SUTRAIce_Exp, SUTRAIce_Power
include("sfcc.jl")

include("sfccsolvers/solvers.jl")

# Extra utilities
"""
    ustrip(x)

Reconstructs the type or function `x` with all numerical quantities stripped of units.
"""
ustrip(x::Union{SFCC,SWRC}) = Flatten.reconstruct(x, map(ustrip, Flatten.flatten(x, Flatten.flattenable, Number)), Number)

include("plotting.jl")

end
