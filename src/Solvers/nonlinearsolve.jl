"""
    SFCCNonlinearSolver{TSolver,TOpts} <: SFCCSolver

SFCCSolver wrapper type for non-linear solvers provided by `NonlinearSolve.jl`. `opts` are keyword
arguments which will be passed to `solve`.
"""
struct SFCCNonlinearSolver{TSolver,TOpts} <: SFCCSolver
    nlsolver::TSolver
    opts::TOpts
    SFCCNonlinearSolver(nlsolver, opts) = new{typeof(nlsolver),typeof(opts)}(nlsolver, opts)
    SFCCNonlinearSolver(nlsolver::TSolver=NewtonRaphson(); tol=1e-3, opts...) where {TSolver} = let opts = tuple(tol, opts...); new{TSolver,typeof(opts)}(nlsolver, opts) end
end

"""
    sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNonlinearSolver, T₀::Number, ::Val{return_all}=Val{true}()) where {return_all}

Solves the given temperature residual objective with `solver` and initial temperature guess `T₀`.
"""
function sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNonlinearSolver, T₀::Number, ::Val{return_all}=Val{true}()) where {return_all}
    resid(T) = obj(T) # extract residual from objective function return value
    f(T, p) = resid.(T)
    u0 = @SVector[T₀]
    prob = NonlinearProblem{false}(f, u0)
    T = first(solve(prob, solver.nlsolver, solver.opts...))
    θw, ∂θw∂T = ∇(T -> obj.f(T; obj.f_kwargs...), T)
    C = obj.hc(θw)
    return return_all ? (; T, θw, C, ∂θw∂T) : T
end
function sfccsolve(obj::SFCCInverseEnthalpyObjective, solver::SFCCNonlinearSolver, T₀::NTuple{2}, ::Val{return_all}=Val{true}()) where {return_all}
    resid(T) = obj(T) # extract residual from objective function return value
    f(T, p) = resid(T)
    prob = IntervalNonlinearProblem{false}(f, T₀)
    T = first(solve(prob, solver.nlsolver, solver.opts...))
    θw, ∂θw∂T = ∇(T -> obj.f(T; obj.f_kwargs...), T)
    C = obj.hc(θw)
    return return_all ? (; T, θw, C, ∂θw∂T) : T
end
