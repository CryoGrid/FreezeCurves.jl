"""
    heatcapacity(c::Number)

Heat capacity function `f(θw,θtot,θsat) = c` where `c` is a pre-specified constant.
"""
heatcapacity(c::Number) = (θw,θtot,θsat) -> c
"""
    heatcapacity(c_frozen::Number, c_thawed::Number; θtot=1.0)

Piecewise constant heat capacity function:
```math
f(θw,θtot,θsat) =
    \\begin{cases}
    c_f & θw < θtot \\\\
    c_t & \\text{otherwise}
    \\end{cases}
```
"""
heatcapacity(c_frozen::Number, c_thawed::Number) = (θw,θtot,θsat) -> IfElse.ifelse(θw < θtot, c_frozen, c_thawed)
