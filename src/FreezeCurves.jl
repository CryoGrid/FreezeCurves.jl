module FreezeCurves

using ForwardDiff
using IfElse
using Interpolations
using IntervalSets
using RecipesBase
using Reexport
using Requires
@reexport using Setfield: @set, @set!
@reexport using Unitful

import Flatten
import Unitful: ustrip

function __init__()
    @require Turing="fce5fe82-541a-59a6-adf8-730c64b5f9a0" begin
        using .Turing
        include("Inference/Inference.jl")
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

include("math.jl")
export SWRC, SoilWaterVolume, BrooksCorey, VanGenuchten
export inflectionpoint
include("swrc.jl")
export SFCC, SFCCSolver, SoilFreezeThawProperties
export PainterKarra, DallAmico, DallAmicoSalt, McKenzie, Westermann, Hu2020, PowerLaw
export SUTRAIce_Exp, SUTRAIce_Power
include("sfcc.jl")
include("Solvers/Solvers.jl")
using .Solvers

# Free water freeze curve
"""
    FreeWater <: FreezeCurve

"Free water" freeze curve in terms of enthalpy (H), total water content (θtot), and
the latent heat of fusion of water (L).
"""
struct FreeWater <: FreezeCurve end
function freewater(H, θtot, L)
    θtot = max(1e-8, θtot)
    Lθ = L*θtot
    I_t = H > Lθ
    I_f = H <= 0.0
    I_c = (H > 0.0) && (H <= Lθ)
    # compute liquid water content -> heat capacity -> temperature
    θw = (I_c*(H/Lθ) + I_t)θtot
    return (;θw, I_t, I_f, I_c, Lθ)
end
(freeW::FreeWater)(H, θtot, L) = freewater(H, θtot, L).θw

# Extra utilities
"""
    ustrip(x)

Reconstructs the type or function `x` with all numerical quantities stripped of units.
"""
ustrip(x::Union{SFCC,SWRC}) = Flatten.reconstruct(x, map(ustrip, Flatten.flatten(x, Flatten.flattenable, Number)), Number)

include("plotting.jl")

end
