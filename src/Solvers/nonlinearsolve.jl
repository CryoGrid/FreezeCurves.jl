using NonlinearSolve

"""
    SFCCNonlinearSolver{TSolver<:AbstractNonlinearSolveAlgorithm,TOpts} <: SFCCSolver

SFCCSolver wrapper type for non-linear solvers provided by `NonlinearSolve.jl`. `opts` are keyword
arguments which will be passed to `solve`.
"""
struct SFCCNonlinearSolver{TSolver<:NonlinearSolve.AbstractNonlinearSolveAlgorithm,TOpts} <: SFCCSolver
    nlsolver::TSolver
    opts::TOpts
    SFCCNonlinearSolver(nlsolver::TSolver; tol=1e-3, opts...) where {TSolver} = let opts = tuple(tol, opts...); new{TSolver,typeof(opts)}(nlsolver, opts) end
end

"""
    sfccsolve(obj::SFCCTemperatureObjective, solver::SFCCNonlinearSolver, T₀::Number)

Solves the given temperature residual objective with `solver` and initial temperature guess `T₀`.
"""
function sfccsolve(obj::SFCCTemperatureObjective, solver::SFCCNonlinearSolver, T₀::Number)
    resid(T) = obj(T) # extract residual from objective function return value
    f(T, p) = resid.(T)
    u0 = @SVector[T₀]
    prob = NonlinearProblem{false}(f, u0)
    return first(solve(prob, solver.nlsolver, solver.opts...))
end
function sfccsolve(obj::SFCCTemperatureObjective, solver::SFCCNonlinearSolver, T₀::NTuple{2})
    resid(T) = obj(T) # extract residual from objective function return value
    f(T, p) = resid(T)
    prob = NonlinearProblem{false}(f, T₀)
    return first(solve(prob, solver.nlsolver, solver.opts...))
end
