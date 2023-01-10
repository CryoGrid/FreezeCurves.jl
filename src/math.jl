"""
    heaviside(x)

Differentiable implementation of heaviside step function, i.e:

``h(x) = \\begin{cases} 1 & x ≥ 0 \\\\ 0 & x < 0 \\end{cases}``
"""
heaviside(x) = IfElse.ifelse(x >= zero(x), 1.0, 0.0)
# Automatic differentiation
∇(f, x) = ∇(typeof(x), f, x)
function ∇(::Type{T}, f, x) where {T}
    res = ForwardDiff.derivative!(ForwardDiff.DiffResult(zero(T), x), f, x)
    return res.value, res.derivs[1]
end
# Function tabulation
"""
    Tabulated(f; kwargs...)

Alias for `tabulate` intended for function types.
"""
Tabulated(f, argknots...) = tabulate(f, argknots...)
"""
    tabulate(f; kwargs...)

Tabulates the given function `f` using a linear, multi-dimensional interpolant.
Knots should be given as keyword arguments `arg = A` where `A` is a `StepRange`
or `Vector` of input values (knots) at which to evaluate the function. `A` may
also be a `Number`, in which case a pseudo-point interpolant will be used
(i.e valid on `[A,A+ϵ]`). Note that all arguments to `f`
must be accepted as keyword arguments.
"""
function tabulate(f; kwargs...)
    initknots(a::AbstractArray) = Interpolations.deduplicate_knots!(a)
    initknots(x::Number) = initknots([x,x])
    interp(::AbstractArray) = Gridded(Linear())
    interp(::Number) = Gridded(Constant())
    extrap(::AbstractArray) = Flat()
    extrap(::Number) = Throw()
    kwargs = (;kwargs...)
    names = keys(kwargs)
    # get knots for each argument, duplicating if only one value is provided
    knots = map(initknots, values(kwargs))
    griddims = map(length, Tuple(knots))
    f_values = zeros(griddims)
    Threads.@threads for inds in collect(Iterators.product(map(N -> 1:N, griddims)...))
        args = map(getindex, knots, inds)
        f_values[inds...] = f(; map(Pair, names, args)...)
    end
    # evaluate function construct interpolant
    f = extrapolate(interpolate(Tuple(knots), f_values, map(interp, values(kwargs))), map(extrap, values(kwargs)))
    return f
end
