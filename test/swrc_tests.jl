using FreezeCurves
using Test
using Unitful

@testset "SWRC" begin
    @testset "VanGenuchten" begin
        f = VanGenuchten()
        let θsat = 0.8,
            θres = 0.0,
            α = 4.0u"m^-1",
            n = 2.0;
            @test isapprox(f(-1e6u"m"; θsat, θres, α, n), 0.0, atol=1e-6)
            @test isapprox(f(0.0u"m"; θsat, θres, α, n), θsat, atol=1e-6)
            θw = f(-0.1u"m"; θsat, θres, α, n)
            @test θw > 0.0 && θw < θsat
            # check inverse
            @test isapprox(Base.inv(f)(θw; θsat, θres, α, n), -0.1u"m", atol=1e-6u"m")
            @test isapprox(Base.inv(f)(θsat; θsat, θres, α, n), 0.0u"m", atol=1e-6u"m")
        end
    end
    @testset "VanGenuchten (units)" begin
        f = ustrip(VanGenuchten(α=0.02u"cm^-1", n=1.4))
        @test f.α == 2.0
        @test f.n == 1.4
    end
    @testset "BrooksCorey" begin
        f = BrooksCorey()
        let θsat = 0.8,
            θres = 0.0,
            λ = 1.0,
            ψₛ = 0.1u"m";
            @test isapprox(f(-1e6u"m"; θsat, θres, ψₛ, λ), 0.0, atol=1e-6)
            @test isapprox(f(-0.1u"m"; θsat, θres, ψₛ, λ), θsat, atol=1e-6)
            @test isapprox(f(-0.01u"m"; θsat, θres, ψₛ, λ), θsat, atol=1e-6)
            θw = f(-0.5u"m"; θsat, θres, ψₛ, λ)
            @test θw > 0.0 && θw < θsat
            # check inverse
            @test isapprox(Base.inv(f)(θw; θsat, θres, ψₛ, λ), -0.5u"m", atol=1e-6u"m")
            @test isapprox(Base.inv(f)(θsat; θsat, θres, ψₛ, λ), -ψₛ, atol=1e-6u"m")
        end
    end
end
