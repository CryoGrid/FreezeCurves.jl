# Automatic differentiation
∇(f, x) = ∇(typeof(x), f, x)
function ∇(::Type{T}, f, x) where {T}
    res = ForwardDiff.derivative!(ForwardDiff.DiffResult(zero(T), x), f, x)
    return res.value, res.derivs[1]
end
"""
    dual(x::Number, ::Type{tag}) where {tag}
    dual(x::A, ::Type{tag}) where {N,T,A<:SVector{N,T},tag}

Constructs a `ForwardDiff.Dual` number (or static array thereof) with `tag` from `x`.
"""
function dual(x::Number, ::Type{tag}) where {tag}
    dualx = ForwardDiff.Dual{tag}(x, one(x))
    return dualx
end
@generated function dual(x::A, ::Type{tag}) where {N,T,A<:SVector{N,T},tag}
    # convert each value of `x` to a ForwardDiff.Dual using `single_seed` to produce the appropriate
    # partial derivatives for each index.
    dual_constructors = (:(ForwardDiff.Dual{tag}(x[$i], ForwardDiff.single_seed(ForwardDiff.Partials{N,eltype(x)}, Val{$i}()))) for i in 1:N)
    return :(SVector{$N}(tuple($(dual_constructors...))))
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
    initknots(a::AbstractArray) = length(a) > 1 ? Interpolations.deduplicate_knots!(a) : Interpolations.deduplicate_knots!(repeat(a,2))
    initknots(x::Number) = initknots([x,x])
    interp(::AbstractArray) = Gridded(Interpolations.Linear())
    interp(::Number) = Gridded(Interpolations.Constant())
    extrap(::AbstractArray) = Interpolations.Flat()
    extrap(::Number) = Interpolations.Throw()
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
