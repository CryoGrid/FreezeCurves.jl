using Enzyme

@testset "Enzyme" begin
    @testset "$name" for (name, (f, x₀)) in autodiff_cases()
        ∂f∂x_reverse = only(Enzyme.gradient(Reverse, f, x₀))
        ∂f∂x_forward = only(Enzyme.gradient(Forward, f, x₀))
        @test isfinite(∂f∂x_reverse)
        @test ∂f∂x_reverse ≈ finitediff(f, x₀) rtol=1e-5
        @test ∂f∂x_forward ≈ ∂f∂x_reverse
    end
end
