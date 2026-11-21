# convenience constants for temperature unit/quantity types
const TemperatureUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TemperatureQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TemperatureUnit}

"""
    enthalpy(T, C, L, θw)

Computes the volumetric enthalpy [J/m^3] given T, C, L, θw.
"""
enthalpy(T, C, L, θw) = T*C + L*θw

"""
    normalize_temperature(x)
    normalize_temperature(x::TemperatureQuantity)

Converts temperature `x` to Kelvin. If `x` has units, `uconvert` is used. Otherwise, if `x` a general numeric type, does a heuristic check `x <= 200`
and assumes the given temperature is in celsius if the condition is satisfied.
"""
function normalize_temperature(x)
    FT = typeof(float(x))
    return x > FT(200) ? x : x + FT(273.15)
end
normalize_temperature(x::TemperatureQuantity) = uconvert(u"K", x)

"""
    heaviside(x)

Differentiable implementation of heaviside step function, i.e:

``h(x) = \\begin{cases} 1 & x ≥ 0 \\\\ 0 & x < 0 \\end{cases}``
"""
heaviside(x) = ifelse(x >= zero(x), one(x), zero(one(x)))

# strip units and conver
stripunits(x::Number) = x
stripunits(x::Quantity) = ustrip(upreferred(x))
stripunits(x::TemperatureQuantity) = ustrip(x)

# convert the number type of x to NF while preserving units
convert_numtype(::Type{NF}, x::Number) where {NF<:Number} = convert(NF, x)
convert_numtype(::Type{NF}, x::Quantity) where {NF<:Number} = convert(NF, ustrip(x))*unit(x)
