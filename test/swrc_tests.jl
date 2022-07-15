using FreezeCurves
using Test
using Unitful

@testset "SWRC" begin
    @testset "VanGenuchten" begin
        f = VanGenuchten()
        let θsat = 0.8,
            α = 4.0u"1/m",
            n = 2.0,
            θres = 0.0;
            @test isapprox(f(-1e6u"m",θsat,θres,α,n), 0.0, atol=1e-6)
            @test isapprox(f(0.0u"m",θsat,θres,α,n), θsat, atol=1e-6)
            θw = f(-0.1u"m",θsat,θres,α,n)
            @test θw > 0.0 && θw < θsat
            # check inverse
            @test isapprox(Base.inv(f)(θw,θsat,θres,α,n), -0.1u"m", atol=1e-6u"m")
        end
    end
end
