using FreezeCurves
using Test
using Unitful

@testset "SFCC" begin
    @testset "McKenzie" begin
        f = McKenzie()
        let θtot = θsat = 0.8,
            γ = 1.0u"K",
            θres = 0.0,
            Tₘ = 0.0u"°C";
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,γ), 0.0, atol=1e-6)
            @test f(0.0u"°C"; θsat,θres,Tₘ,γ) ≈ θtot
            θw = f(-0.1u"°C"; θsat,θres,Tₘ,γ)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "Westermann freeze curve" begin
        f = Westermann()
        let θtot = θsat = 0.8,
            δ = 0.1u"K",
            θres = 0.0,
            Tₘ = 0.0u"°C";
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,δ), 0.0, atol=1e-2)
            @test f(0.0u"°C"; θsat,θres,Tₘ,δ) ≈ θtot
            θw = f(-0.1u"°C"; θsat,θres,Tₘ,δ)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "Hu2020 freeze curve" begin
        f = Hu2020()
        let θtot = θsat = 0.8,
            b = 0.001,
            θres = 0.0,
            Tₘ = 0.0u"°C";
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,b), 0.0, atol=1e-2)
            @test f(0.0u"°C"; θsat,θres,Tₘ,b) ≈ θtot
            θw = f(-0.1u"°C"; θsat,θres,Tₘ,b)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "PowerLaw freeze curve" begin
        f = PowerLaw()
        let θtot = θsat = 0.8,
            α = 0.01,
            β = 0.5,
            θres = 0.0;
            @test isapprox(f(-10.0u"°C"; θsat,θres,α,β), 0.0, atol=1e-2)
            @test f(0.0u"°C"; θsat,θres,α,β) ≈ θtot
            θw = f(-0.5u"°C"; θsat,θres,α,β)
            @test θw > 0.0 && θw < θtot
        end
    end
    @testset "PainterKarra freeze curve" begin
        f = PainterKarra(β=1.0, ω=0.1)
        let θtot = θsat = 0.8,
            α = 4.0u"1/m",
            n = 2.0,
            Tₘ = 0.0u"°C",
            θres = 0.0;
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,α,n), 0.0, atol=1e-3)
            @test f(0.0u"°C"; θsat,θres,Tₘ,α,n) ≈ θtot
            θw = f(-0.1u"°C"; θsat,θres,Tₘ,α,n)
            @test θw > 0.0 && θw < θtot
            res = f(-0.1u"°C", 1.0, Val{:all}(); θsat,θres,Tₘ,α,n)
            @test keys(res) == (:θw,:ψ,:Tstar)
            ψ = f(-0.1u"°C", 1.0, Val{:ψ}(); θsat,θres,Tₘ,α,n)
            @test ψ == res.ψ
            Tstar = f(-0.1u"°C", 1.0, Val{:Tstar}(); θsat,θres,Tₘ,α,n)
            @test Tstar == res.Tstar
        end
    end
    @testset "DallAmico freeze curve" begin
        f = DallAmico()
        let θtot = θsat = 0.8,
            α = 4.0u"1/m",
            n = 2.0,
            Tₘ = 0.0u"°C",
            θres = 0.0;
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,α,n), 0.0, atol=1e-3)
            @test f(0.0u"°C"; θsat,θres,Tₘ,α,n) ≈ θtot
            θw = f(-0.1u"°C"; θsat,θres,Tₘ,α,n)
            @test θw > 0.0 && θw < θtot
            pkfc = PainterKarra(swrc=VanGenuchten(α=α, n=n))
            @test f(-0.1u"°C"; θsat,θres,Tₘ,α,n) ≈ pkfc(-0.1u"°C"; θsat,θres,Tₘ,α,n)
        end
    end
    @testset "DallAmicoSalt freeze curve" begin
        f = DallAmicoSalt()
        let θtot = θsat = 0.8,
            α = 4.0u"1/m",
            n = 2.0,
            Tₘ = 0.0u"°C",
            saltconc = 800.0u"mol/m^3",
            θres = 0.0;
            @test isapprox(f(-10.0u"°C"; θsat,θres,Tₘ,saltconc,α,n), 0.0, atol=1e-3)
            @test f(-0.1u"°C"; θsat,θres,Tₘ,saltconc,α,n) ≈ θtot
            θw = f(-5.0u"°C"; θsat,θres,Tₘ,saltconc,α,n)
            @test θw > 0.0 && θw < θtot
            θw_nosalt = f(-5.0u"°C"; θsat,θres,Tₘ,saltconc=zero(saltconc),α,n)
            @test θw > θw_nosalt
            res = f(-0.1u"°C", 1.0, Val{:all}(); θsat,θres,Tₘ,α,n)
            @test keys(res) == (:θw,:ψ,:Tstar)
            ψ = f(-0.1u"°C", 1.0, Val{:ψ}(); θsat,θres,Tₘ,α,n)
            @test ψ == res.ψ
            Tstar = f(-0.1u"°C", 1.0, Val{:Tstar}(); θsat,θres,Tₘ,α,n)
            @test Tstar == res.Tstar
        end
    end
end
