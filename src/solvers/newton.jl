"""
    SFCCNewtonSolver <: SFCCSolver

Specialized implementation of Newton's method with backtracking line search for resolving
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
function sfccsolve(solver::SFCCNewtonSolver, f::F, f_args, f_hc, L, H, T₀::Nothing=nothing) where {F}
    T₀ = IfElse.ifelse(H < zero(H), H / f_hc(0.0), zero(H))
    return sfccsolve(solver, f, f_args, f_hc, L, H, T₀)
end
using ForwardDiff
function sfccsolve(solver::SFCCNewtonSolver, f::F, f_args, f_hc, L, H, T₀) where {F}
    fc(T) = f(T, f_args...)
    resid(T) = temperature_residual(f, f_args, f_hc, L, H, T)
    T = T₀
    α₀ = solver.α₀
    τ = solver.τ
    # compute initial residual
    Tres, θw, C = resid(T)
    itercount = 0
    T_converged = false
    while abs(Tres) > solver.abstol && abs(Tres) / abs(T) > solver.reltol
        if itercount >= solver.maxiter
            return (;T, Tres, θw, itercount, T_converged)
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
                @warn "Backtracking failed; this should not happen. Current state: α=$α, T=$T, T̂=$T̂, residual $(T̂res), initial residual: $(Tres)"
                return (;T, Tres, θw, itercount, T_converged)
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
    T_converged = true
    return (;T, Tres, θw, itercount, T_converged)
end
