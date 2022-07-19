module Solvers

    using ..FreezeCurves
    using ..FreezeCurves: ∇

    using ForwardDiff
    using IfElse
    using Interpolations
    using NLsolve
    using StaticArrays
    using Requires
    using Unitful

    export sfccsolve, SFCCTemperatureObjective

    """
        AbstractSFCCObjective

    Base type for non-linear SFCC optimization objectives.
    """
    abstract type AbstractSFCCObjective end
    """
        SFCCTemperatureObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH} <: AbstractSFCCObjective

    Optimization objective for finding a temperature, `T`, wich resolves the conservation equation:
    ```math
    0 = (H - Lθ)/C - T
    ```
    given a fixed enthalpy value `H` and freeze curve arguments `f_args`.
    """
    Base.@kwdef struct SFCCTemperatureObjective{TF,Tkwargs<:NamedTuple,Thc,TL,TH} <: AbstractSFCCObjective
        f::TF
        f_kwargs::Tkwargs
        hc::Thc
        L::TL
        H::TH
    end
    @inline (obj::SFCCTemperatureObjective)(T) = FreezeCurves.temperature_residual(obj.f, obj.f_args, obj.hc, obj.L, obj.H, T, true)

    """
        sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀)

    Solve the given objective `obj` using `solver` and initial guess `x₀`.
    """
    sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀) = error("not implemented for objective $(typeof(obj)) with solver $(typeof(solver))")

    export heatcapacity
    include("heatcap.jl")
    export SFCCNewtonSolver
    include("newton.jl")
    export SFCCPreSolver
    include("presolver.jl")
    # require NonlinearSolve.jl for generic nonlinear solver
    @require NonlinearSolve="8913a72c-1f9b-4ce2-8d82-65094dcecaec" begin
        export SFCCNonlinearSolver
        include("nonlinearsolve.jl")
    end
end
