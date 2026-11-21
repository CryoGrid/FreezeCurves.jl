using ForwardDiff

@testset "ForwardDiff" begin
    @testset "$name" for (name, (f, x₀)) in autodiff_cases()
        ∂f∂x = ForwardDiff.derivative(f, x₀)
        @test isfinite(∂f∂x)
        @test ∂f∂x ≈ finitediff(f, x₀) rtol=1e-5
    end
end
