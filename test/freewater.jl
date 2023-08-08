using FreezeCurves
using Test

@testset "Free water freeze curve" begin
    fc = FreeWater()
    L = 3.34e8u"J/m^3"
    θtot = 0.5
    θres = 0.0
    @test fc(-1e6u"J/m^3"; θtot, θres, L) ≈ θres
    @test fc(0.0u"J/m^3"; θtot, θres, L) ≈ θres
    @test fc(L*θtot/2; θtot, θres, L) ≈ θtot/2
    @test fc(L*θtot; θtot, θres, L) ≈ θtot
    @test fc(2*L*θtot; θtot, θres, L) ≈ θtot
end
