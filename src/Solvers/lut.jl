"""
    LUT{N,T,coordnames,TCoords,TData}

Generic implementation of an N-dimensional lookup table with support for reverse lookups
w.r.t one variable.
"""
struct LUT{N,T,coordnames,TCoords,TData<:AbstractArray{T,N}}
    coords::NamedTuple{coordnames,TCoords}
    interp::TData
end
function (lut::LUT)(; coords...)
    @assert all(map(∈(keys(coords)), keys(lut.coords))) "All coordinates must be specified for forward lookup."
    queries = map(k -> coords[k], keys(lut.coords))
    return lut.interp(queries...)
end
@generated function free_coord_index(lut::LUT{N,T,coordnames}; kwqueries...) where {N,T,coordnames}
    querynames(::Type{<:Iterators.Pairs{TK,TV,TT,<:NamedTuple{names}}}) where {TK,TV,TT,names} = names
    i = findfirst(∉(querynames(kwqueries)), coordnames)
    return :($i)
end
function coordinatesearch(f, A::AbstractArray, query)
    # prob = IntervalNonlinearProblem{false}((u,p) -> f(u) - query, extrema(A))
    # ustar = solve(prob, Falsi()).u
    # res = searchsorted(A, ustar)
    # hi, lo = first(res), last(res)
    # f_lo = f(A[lo])
    # f_hi = f(A[hi])
    # return A[lo] + (query - f_lo)*(A[hi] - A[lo])/(f_hi - f_lo)
    lo = 1
    hi = length(A)
    f_lo = f(A[lo])
    f_hi = f(A[hi])
    if query < f_lo
        return f_lo + query*(f_lo - f(A[lo+1]))/(A[lo] - A[lo+1])
    elseif query > f_hi
        return f_hi + query*(f_hi - f(A[hi-1]))/(A[hi] - A[hi-1])
    else
        while hi - lo > 1
            mid = lo + Int(floor((hi - lo) / 2))
            f_mid = f(A[mid])
            if query <= f_mid
                hi = mid
            else
                lo = mid
            end
        end
        f_lo = f(A[lo])
        f_hi = f(A[hi])
        return A[lo] + (query - f_lo)*(A[hi] - A[lo])/(f_hi - f_lo)
    end
end
function (lut::LUT{N,T})(value; kwqueries...) where {N,T}
    @assert length(kwqueries) == N - 1 "Inverse lookup in $(N)D LUT requires $(N-1) coordinates to be specified."
    target_coord_idx = free_coord_index(lut; kwqueries...)
    queries = map(keys(lut.coords)) do k
        get(kwqueries, k, :)
    end |> NamedTuple{keys(lut.coords)}
    f(x) = lut.interp(replace(values(queries), Colon() => x)...)
    return coordinatesearch(f, lut.coords[target_coord_idx], value)
end
"""
    build_lut(f, coords::Pair{Symbol,<:Union{Number,AbstractVector}}...; f_kwargs...)

Builds an n-dimensional lookup table `LUT` by invoking `f` on the given `coords`. `f` should
be a function which accepts keyword arguments matching `coords` as well as any additional
constant keyword arguments `f_kwargs`.

Example:
```julia
f(; x=1.0, y=1.0) = x + y
lut = build_lut(f, :x => -10.0:10.0, :y=-10.0:10.0)
lut(; x=2.2, y=1.5)
# output
3.7
```
"""
function build_lut(f, coords::Pair{Symbol,<:Union{Number,AbstractVector}}...; f_kwargs...)
    @assert all(map(issorted ∘ last, coords))
    lut_f(; kwcoords...) = f(; kwcoords..., f_kwargs...)
    # make sure every coordinate range is an array
    new_coords = map(kv -> first(kv) => vec(collect(last(kv))), coords)
    interp = FreezeCurves.tabulate(lut_f; new_coords...)
    return LUT((; new_coords...), interp)
end