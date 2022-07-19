"""
    heatcapacity(c::Number)

Heat capacity function `f(H,T,θw) = c` where `c` is a pre-specified constant.
"""
heatcapacity(c::Number) = θw -> c
"""
    heatcapacity(c_frozen::Number, c_thawed::Number; θtot=1.0)

Piecewise constant heat capacity function:
```math
f(θw) =
    \\begin{cases}
    c_f & θw < θtot \\\\
    c_t & \\text{otherwise}
    \\end{cases}
```
"""
heatcapacity(c_frozen::Number, c_thawed::Number; θtot=1.0) = θw -> IfElse.ifelse(θw < θtot, c_frozen, c_thawed)
