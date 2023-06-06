using FreezeCurves
using FreezeCurves: SFCCPreSolverCacheND, SFCCPreSolverCache1D
using NonlinearSolve
using Test

using Base.Iterators

@testset "temperature_residual" begin
    freezecurves = (ustrip(McKenzie()), ustrip(Westermann()), ustrip(DallAmico()))
    L = 3.34e8
    T₀ = -1.0
    for fc in freezecurves
        vol = SoilWaterVolume(fc)
        θsat = vol.θsat
        hc = FreezeCurves.heatcapacity(1.9e6, 4.2e6)
        θw_true = fc(-0.01; θsat)
        C_true = hc(θw_true, θsat, θsat)
        H_true = FreezeCurves.enthalpy(-0.01, C_true, L, θw_true)
        T_resid, θw, C = @inferred FreezeCurves.temperature_residual(fc, hc, L, H_true, -0.1, 1.0)
        @test abs(T_resid) > 0
    end
end

@testset "Solvers" begin
    freezecurves = (ustrip(McKenzie()), ustrip(Westermann()), ustrip(DallAmico()), ustrip(PainterKarra()))
    T_cases = [1.0,-0.01,-0.1,-10.0]
    L = 3.34e8
    T₀ = -1.0
    @testset "Newton" begin
        solver = SFCCNewtonSolver(abstol=1e-5)
        for (T_true, fc) in product(T_cases, freezecurves)
            vol = SoilWaterVolume(fc)
            θsat = vol.θsat
            hc = FreezeCurves.heatcapacity(1.9e6, 4.2e6)
            θw_true = fc(T_true; θsat)
            C_true = hc(θw_true, θsat, θsat)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCInverseEnthalpyObjective(fc, (; θsat), hc, L, H_true, 1.0)
            res = sfccsolve(obj, solver, T₀)
            @test isapprox(res.T, T_true, atol=1e-3)
        end
    end
    @testset "Nonlinear" begin
        # we use Falsi bracketing method for testing here because NewtonRaphson seems to have
        # convergence issues for some freeze curves and initial guesses.
        solver = SFCCNonlinearSolver(Falsi(), abstol=1e-5)
        for (T_true, fc) in product(T_cases, freezecurves)
            vol = SoilWaterVolume(fc)
            θsat = vol.θsat
            hc = FreezeCurves.heatcapacity(1.9e6, 4.2e6)
            θw_true = fc(T_true; θsat)
            C_true = hc(θw_true, θsat, θsat)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCInverseEnthalpyObjective(fc, (; θsat), hc, L, H_true, 1.0)
            # bracketing method requires a tuple range as initial guess
            res = sfccsolve(obj, solver, (-20.0,10.0))
            @test isapprox(res.T, T_true, atol=1e-3)
        end
    end
    @testset "Presolver" begin
        solver = SFCCPreSolver(SFCCPreSolverCache1D(), errtol=1e-5)
        for (T_true, fc) in product(T_cases, freezecurves)
            vol = SoilWaterVolume(fc)
            θsat = vol.θsat
            hc = FreezeCurves.heatcapacity(1.9e6, 4.2e6)
            FreezeCurves.initialize!(solver, fc, hc; θsat)
            θw_true = fc(T_true; θsat)
            C_true = hc(θw_true, θsat, θsat)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCInverseEnthalpyObjective(fc, (; θsat), hc, L, H_true, 1.0)
            res = sfccsolve(obj, solver, T₀)
            @test isapprox(res.T, T_true, atol=1e-3)
        end
        solver = SFCCPreSolver(SFCCPreSolverCacheND())
        for (T_true, fc) in product(T_cases, freezecurves)
            vol = SoilWaterVolume(fc)
            θsat = vol.θsat
            hc = FreezeCurves.heatcapacity(1.9e6, 4.2e6)
            FreezeCurves.initialize!(solver, fc, hc; θsat)
            sat = 0.95
            θw_true = fc(T_true, sat; θsat)
            C_true = hc(θw_true, θsat*sat, θsat)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCInverseEnthalpyObjective(fc, (; θsat), hc, L, H_true, sat)
            res = sfccsolve(obj, solver, T₀)
            @test isapprox(res.T, T_true, atol=0.1)
        end
    end
end
