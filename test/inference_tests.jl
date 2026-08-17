using FreezeCurves
using Turing
using Test

@testset "SFCC inference" begin

    Trange = vcat(-5.0u"°C":0.1u"K":-0.11u"°C", -0.1u"°C":0.001u"K":0.0u"°C")

    @testset "Dall'Amico" begin
        fc = DallAmico(swrc=VanGenuchten(α=0.1u"1/m", n=1.8))
        θtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))
        sfcc_model = SFCCModel(fc)
        # condition on data
        m = sfcc_model(Trange, θtrue)
        _, pred = m()
        @test all(isfinite.(pred.θw))
        @test all(0 .< pred.θw .< 1)
    end

    @testset "PainterKarra" begin
        fc = PainterKarra(swrc=VanGenuchten(α=0.1u"1/m", n=1.8))
        θtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))
        sfcc_model = SFCCModel(fc)
        # condition on data
        m = sfcc_model(Trange, θtrue)
        _, pred = m()
        @test all(isfinite.(pred.θw))
        @test all(0 .< pred.θw .< 1)
    end

    @testset "Hu2020" begin
        fc = Hu2020()
        θtrue = min.(0.5, max.(fc.(Trange) .+ randn(length(Trange)).*0.02, 0.0))
        sfcc_model = SFCCModel(fc)
        # condition on data
        m = sfcc_model(Trange, θtrue)
        _, pred = m()
        @test all(isfinite.(pred.θw))
        @test all(0 .< pred.θw .< 1)
    end
end
