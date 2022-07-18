module Solvers

    using ..FreezeCurves
    using StaticArrays
    using Unitful

    export sfccsolve, SFCCTemperatureObjective

    """
        AbstractSFCCObjective

    Base type for non-linear SFCC optimization objectives.
    """
    abstract type AbstractSFCCObjective end
    """
        SFCCTemperatureObjective{TF,Targs,Thc,TL,TH} <: AbstractSFCCObjective

    Optimization objective for finding a temperature, `T`, wich resolves the conservation equation:
    ```math
    0 = (H - Lθ)/C - T
    ```
    given a fixed enthalpy value `H` and freeze curve arguments `f_args`.
    """
    struct SFCCTemperatureObjective{TF,Targs,Thc,TL,TH} <: AbstractSFCCObjective
        f::TF
        f_args::Targs
        f_hc::Thc
        L::TL
        H::TH
    end
    @inline (obj::SFCCTemperatureObjective)(T) = FreezeCurves.temperature_residual(obj.f, obj.f_args, obj.f_hc, obj.L, obj.H, T, true)

    """
        sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀)

    Solve the given objective `obj` using `solver` and initial guess `x₀`.
    """
    sfccsolve(obj::AbstractSFCCObjective, solver::SFCCSolver, x₀) = error("not implemented for objective $(typeof(obj)) with solver $(typeof(solver))")

    export heatcapacity
    include("heatcap.jl")
    export SFCCNewtonSolver
    include("newton.jl")
    export SFCCNonlinearSolver
    include("nlsolve.jl")
    export SFCCPreSolver
    include("presolver.jl")
end
