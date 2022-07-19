module FreezeCurves

using ForwardDiff
using IfElse
using RecipesBase
using Reexport
using Requires
@reexport using Setfield: @set, @set!
@reexport using Unitful

import Flatten
import Unitful: ustrip

const TemperatureUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TemperatureQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TemperatureUnit}

"""
    FreezeCurve

Base type for freezing characteristic curves which relate temperature or enthalpy to unfrozen water content.
"""
abstract type FreezeCurve end

# Free water freeze curve
"""
    FreeWater <: FreezeCurve

"Free water" freeze curve in terms of enthalpy (H), total water content (θwi), and
the latent heat of fusion of water (L).
"""
struct FreeWater <: FreezeCurve end
function freewater(H, L, θwi)
    θwi = max(1e-8, θwi)
    Lθ = L*θwi
    I_t = H > Lθ
    I_f = H <= 0.0
    I_c = (H > 0.0) && (H <= Lθ)
    # compute liquid water content -> heat capacity -> temperature
    θw = (I_c*(H/Lθ) + I_t)θwi
    return (;θw, I_t, I_f, I_c, Lθ)
end
(freeW::FreeWater)(H, θwi, L) = freewater(H, θwi, L).θw

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
export SWRC, SoilWaterProperties, VanGenuchten
include("swrc.jl")
export SFCC, SFCCFunction, SFCCSolver, SoilFreezeThawProperties, DallAmico, DallAmicoSalt, McKenzie, Westermann
include("sfcc.jl")
include("Solvers/Solvers.jl")

# Extra utilities
"""
    ustrip(x)

Reconstructs the type or function `x` with all numerical quantities stripped of units.
"""
ustrip(x::Union{SFCCFunction,SWRCFunction}) = Flatten.reconstruct(x, map(ustrip, Flatten.flatten(x, Flatten.flattenable, Number)), Number)

include("plotting.jl")

end
