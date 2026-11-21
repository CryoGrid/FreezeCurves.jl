include("ad_test_cases.jl")

@testset "Automatic differentiation" begin
    include("forwarddiff_tests.jl")
    include("enzyme_tests.jl")
    include("zygote_tests.jl")
end
