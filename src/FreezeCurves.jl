module FreezeCurves

using ForwardDiff
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
        export SFCCModel, sfccpriors, swvpriors
        include("inference/inference.jl")
    end
    # require NonlinearSolve.jl for generic nonlinear solver
    @require NonlinearSolve="8913a72c-1f9b-4ce2-8d82-65094dcecaec" begin
        using .NonlinearSolve
        export SFCCNonlinearSolver
        include("sfccsolvers/nonlinearsolve.jl")
    end
end

include("utils.jl")

include("math.jl")

"""
    FreezeCurve

Base type for freezing characteristic curves which relate temperature or enthalpy to unfrozen water content.
"""
abstract type FreezeCurve end

export FreeWater
include("freewater.jl")

export SWRC, BrooksCorey, VanGenuchten, watercontent
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

Reconstructs the type or function `x` with all numerical quantities stripped of units. All units are first converted
to preferred SI units before being stripped.
"""
ustrip(x::Union{SFCC,SWRC}) = Flatten.reconstruct(x, map(stripunits, Flatten.flatten(x, Flatten.flattenable, Number)), Number)

"""
    convert_numtype(::Type{NF}, x::Union{SFCC,SWRC,SoilFreezeThawProperties}) where {NF<:Number}

Reconstructs the type or function `x` with all numerical quantities (including those in nested structs)
converted to the number type `NF`. Units are preserved.
"""
function convert_numtype(::Type{NF}, x::Union{SFCC,SWRC,SoilFreezeThawProperties}) where {NF<:Number}
    values = map(Base.Fix1(convert_numtype, NF), Flatten.flatten(x, Flatten.flattenable, Number))
    return Flatten.reconstruct(x, values, Number)
end

include("plotting.jl")

end
