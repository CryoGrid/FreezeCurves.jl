using FreezeCurves
using FreezeCurves.Solvers
using Test

using Base.Iterators

@testset "Solvers" begin
    freezecurves = (ustrip(McKenzie()), ustrip(Westermann()), ustrip(DallAmico()))
    T_cases = [1.0,-0.001,-0.1,-10.0]
    L = 3.34e8
    T₀ = -1.0
    @testset "Newton" begin
        solver = SFCCNewtonSolver(abstol=1e-4)
        for (T_true, fc) in product(T_cases, freezecurves)
            props = SoilWaterProperties(fc)
            θtot, θsat = props.θtot, props.θsat
            hc = heatcapacity(1.9e6, 4.2e6; θtot)
            θw_true = fc(T_true; θtot, θsat)
            C_true = hc(θw_true)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCTemperatureObjective(fc, (; θtot, θsat), hc, L, H_true)
            res = sfccsolve(obj, solver, T₀, return_all=true)
            @test isapprox(res.T, T_true, atol=1e-3)
        end
    end
    @testset "Presolver" begin
        solver = SFCCPreSolver()
        for (T_true, fc) in product(T_cases, freezecurves)
            props = SoilWaterProperties(fc)
            θtot, θsat = props.θtot, props.θsat
            hc = heatcapacity(1.9e6, 4.2e6; θtot)
            Solvers.initialize!(solver, fc, hc, θtot, θsat; L)
            θw_true = fc(T_true; θtot, θsat)
            C_true = hc(θw_true)
            H_true = FreezeCurves.enthalpy(T_true, C_true, L, θw_true)
            # use default freeze curve parameter values
            obj = SFCCTemperatureObjective(fc, (; θtot, θsat), hc, L, H_true)
            res = sfccsolve(obj, solver, T₀, return_all=true)
            @test isapprox(res.T, T_true, atol=1e-3)
        end
    end
end
