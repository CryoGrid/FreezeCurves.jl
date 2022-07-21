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
    Tabulated(f, argknots...)

Alias for `tabulate` intended for function types.
"""
Tabulated(f, argknots...) = tabulate(f, argknots...)
"""
    tabulate(f, argknots::Pair{Symbol,<:Union{Number,AbstractArray}}...)

Tabulates the given function `f` using a linear, multi-dimensional interpolant.
Knots should be given as pairs `:arg => A` where `A` is a `StepRange` or `Vector`
of input values (knots) at which to evaluate the function. `A` may also be a
`Number`, in which case a pseudo-point interpolant will be used (i.e valid on
`[A,A+ϵ]`). No extrapolation is provided by default but can be configured via
`Interpolations.extrapolate`.
"""
function tabulate(f, argknots::Pair{Symbol,<:Union{Number,AbstractArray}}...)
    initknots(a::AbstractArray) = Interpolations.deduplicate_knots!(a)
    initknots(x::Number) = initknots([x,x])
    interp(::AbstractArray) = Gridded(Linear())
    interp(::Number) = Gridded(Constant())
    extrap(::AbstractArray) = Flat()
    extrap(::Number) = Throw()
    names = map(first, argknots)
    # get knots for each argument, duplicating if only one value is provided
    knots = map(initknots, map(last, argknots))
    f_argnames = argnames(f)
    @assert all(map(name -> name ∈ names, f_argnames)) "Missing one or more arguments $f_argnames in $f"
    arggrid = Iterators.product(knots...)
    # evaluate function construct interpolant
    f = extrapolate(interpolate(Tuple(knots), map(Base.splat(f), arggrid), map(interp ∘ last, argknots)), map(extrap ∘ last, argknots))
    return f
end
