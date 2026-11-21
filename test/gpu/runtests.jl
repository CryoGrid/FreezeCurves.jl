# GPU tests; these are not part of the default test suite and require a CUDA-capable GPU.
# Run with: julia --project=test/gpu -e 'using Pkg; Pkg.instantiate(); include("test/gpu/runtests.jl")'
using CUDA
using FreezeCurves
using Test

CUDA.functional() || error("CUDA is not functional; GPU tests require a CUDA-capable GPU")
CUDA.allowscalar(false)

const NF = Float32

# evaluate f elementwise over CPU and GPU arrays and check that the results agree and stay in NF
function test_gpu_map(f, args...)
    cpu_result = map(f, args...)
    gpu_result = map(f, map(cu, args)...)
    @test gpu_result isa CuArray{NF}
    @test Array(gpu_result) ≈ cpu_result
end

@testset "GPU" begin
    T = collect(range(NF(-10), NF(1), length=1024))
    sat = collect(range(NF(0.5), NF(1), length=1024))
    θsat = collect(range(NF(0.3), NF(0.6), length=1024))
    @testset "SWRC: $(nameof(typeof(swrc)))" for swrc in (VanGenuchten(NF), BrooksCorey(NF))
        f = ustrip(swrc)
        ψ = collect(range(NF(-100), NF(1), length=1024))
        θ = collect(range(NF(0.05), NF(0.5), length=1024))
        ∂f∂ψ = FreezeCurves.derivative(f)
        f⁻¹ = Base.inv(f)
        test_gpu_map((ψ, θsat) -> f(ψ; θsat), ψ, θsat)
        test_gpu_map((ψ, θsat) -> ∂f∂ψ(ψ; θsat), ψ, θsat)
        test_gpu_map(θ -> f⁻¹(θ; θsat=NF(0.5)), θ)
    end
    freezecurves = (
        PainterKarra(NF),
        DallAmico(NF),
        DallAmicoSalt(NF),
        McKenzie(NF),
        Westermann(NF),
        Langer(NF),
        Hu2020(NF),
        PowerLaw(NF),
    )
    @testset "SFCC: $(nameof(typeof(fc)))" for fc in freezecurves
        f = ustrip(fc)
        test_gpu_map((T, sat, θsat) -> f(T, sat; θsat), T, sat, θsat)
    end
    @testset "SFCC with SWRC parameters: $(nameof(typeof(fc)))" for fc in (DallAmico(NF), PainterKarra(NF))
        f = ustrip(fc)
        test_gpu_map(T -> f(T, one(NF), Val{:ψ}(); watercontent(NF)..., α=NF(4), n=NF(1.5)), T)
        test_gpu_map(T -> f(T, one(NF), Val{:Tstar}(); watercontent(NF)...), T)
    end
end
