"""
    SFCCNewtonSolver <: SFCCSolver

Fast, specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to
discontinuities and strong non-linearity in most common soil freeze curves.
"""
Base.@kwdef struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 100 # maximum number of iterations
    abstol::Float64 = 1e-2 # absolute tolerance for convergence
    reltol::Float64 = 1e-4 # relative tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.7 # step size decay for backtracking
end
# Newton solver implementation
"""
    sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNewtonSolver, T₀::Number, ::Val{return_all}=Val{true}(); error_when_not_converged=true)

Solves `obj` using the specialized Newton `solver` and returns the result. If `return_all` is `true`,
a named tuple `(;T, Tres, θw, C, itercount)` is returned; otherwise (by default), only the temperature solution is returned.
"""
function sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNewtonSolver, T₀::Number, ::Val{return_all}=Val{true}(); error_when_not_converged=true) where {return_all}
    resid(T) = FreezeCurves.temperature_residual(obj.f, obj.f_kwargs, obj.hc, obj.L, adstrip(obj.H), T, adstrip(obj.sat))
    T = adstrip(T₀)
    α₀ = solver.α₀
    τ = solver.τ
    # compute initial residual
    Tres, θw, C = resid(T)
    itercount = 0
    while abs(Tres) > solver.abstol && abs(Tres) / abs(T) > solver.reltol
        if itercount >= solver.maxiter
            iterstate = (;T, θw, C, ∂θw∂T=NaN, Tres, itercount)
            msg = "Failed to converge within $(solver.maxiter) iterations: $iterstate"
            if error_when_not_converged
                error(msg)
            else
                @warn msg
                return iterstate
            end
            return iterstate
        end
        # derivative of residual
        ∂Tres∂T = ForwardDiff.derivative(first ∘ resid, T)
        α = α₀ / ∂Tres∂T
        T̂ = T - α*Tres
        # do first residual check outside of loop;
        # this way, we don't decrease α unless we have to.
        T̂res, θw, C = resid(T̂)
        inneritercount = 0
        # simple backtracking line search to avoid jumping over the solution
        while sign(T̂res) != sign(Tres)
            if inneritercount > 100
                error("Backtracking failed; this should not happen. There is a bug! Current state: α=$α, T=$T, T̂=$T̂, residual $(T̂res), initial residual: $(Tres)")
            end
            α = α*τ # decrease step size by τ
            T̂ = T - α*Tres # new guess for T
            T̂res, θw, C = resid(T̂)
            inneritercount += 1
        end
        T = T̂ # update T
        Tres = T̂res # update residual
        itercount += 1
    end
    # Re-evaluate freeze curve at solution
    θw, ∂θw∂T = ∇(T -> obj.f(T; obj.f_kwargs...), T)
    return return_all ? (;T, θw, C, ∂θw∂T, Tres, itercount) : T
end
