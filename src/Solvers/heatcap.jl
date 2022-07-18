"""
    heatcapacity(f::Function, ::Val{var})

Heat capacity function constructor which creates a function `(H,T,θw) -> f(var)` where
var is one of `:T`, `:H`, or `:θw`. This is a convenience function for cases where the
user has a heat capacity function which depends solely on enthalpy, temperature, or
unfrozen water content.
"""
heatcapacity(f::Function, ::Val{:H}) = (H,T,θw) -> f(H)
heatcapacity(f::Function, ::Val{:T}) = (H,T,θw) -> f(T)
heatcapacity(f::Function, ::Val{:θw}) = (H,T,θw) -> f(θw)
"""
    heatcapacity(c::Number)

Heat capacity function `f(H,T,θw) = c` where `c` is a pre-specified constant.
"""
heatcapacity(c::Number) = (H,T,θw) -> c
"""
    heatcapacity(c_frozen::Number, c_thawed::Number; T_melt::Number = 0.0u"°C")

Piecewise constant heat capacity function:
```math
f(H,T,θw) =
    \\begin{cases}
    c_f & T \\leq T_m \\\\
    c_t & \\text{otherwise}
    \\end{cases}
```
"""
heatcapacity(c_frozen::Number, c_thawed::Number; T_melt::Number = 0.0u"°C") = (H,T,θw) -> IfElse.ifelse(T <= T_melt, c_frozen, c_thawed)
