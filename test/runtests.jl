using FreezeCurves
using Test
using Unitful

# test groups can be selected via test_args, e.g. Pkg.test(test_args=["autodiff"]); all groups run by default
const TEST_GROUPS = isempty(ARGS) ? ["core", "autodiff"] : ARGS

@testset "FreezeCurves.jl" begin
    if "core" in TEST_GROUPS
        include("freewater_tests.jl")
        include("swrc_tests.jl")
        include("sfcc_tests.jl")
        include("solver_tests.jl")
        include("inference_tests.jl")
    end
    if "autodiff" in TEST_GROUPS
        include("autodiff/runtests.jl")
    end
end
