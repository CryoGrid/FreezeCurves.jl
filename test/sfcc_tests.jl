using FreezeCurves
using Test
using Unitful

@testset "SFCC" begin
    @testset "McKenzie" begin
        f = McKenzie()
        let θtot = 0.8,
            γ = 1.0u"K",
            θres = 0.0,
            Tₘ = 0.0u"°C";
            @test isapprox(f(-10.0u"°C",θtot,θtot,θres,Tₘ,γ), 0.0, atol=1e-6)
            @test f(0.0u"°C",θtot,θtot,θres,Tₘ,γ) ≈ θtot
            θw = f(-0.1u"°C",θtot,θtot,θres,Tₘ,γ)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "Westermann freeze curve" begin
        f = Westermann()
        let θtot = 0.8,
            δ = 0.1u"K",
            θres = 0.0,
            Tₘ = 0.0u"°C";
            @test isapprox(f(-10.0u"°C",θtot,θtot,θres,Tₘ,δ), 0.0, atol=1e-2)
            @test f(0.0u"°C",θtot,θtot,θres,Tₘ,δ) ≈ θtot
            θw = f(-0.1u"°C",θtot,θtot,θres,Tₘ,δ)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "DallAmico freeze curve" begin
        f = DallAmico()
        let θtot = 0.8,
            α = 4.0u"1/m",
            n = 2.0,
            Tₘ = 0.0u"°C",
            θres = 0.0;
            @test isapprox(f(-10.0u"°C",θtot,θtot,θres,Tₘ,α,n), 0.0, atol=1e-3)
            @test f(0.0u"°C",θtot,θtot,θres,Tₘ,α,n) ≈ θtot
            θw = f(-0.1u"°C",θtot,θtot,θres,Tₘ,α,n)
            @test θw > 0.0 && θw < θtot
        end
    end
end
