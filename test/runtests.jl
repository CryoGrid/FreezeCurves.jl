using FreezeCurves
using Test
using Unitful

@testset "FreezeCurves.jl" begin
    include("freewater_tests.jl")
    include("swrc_tests.jl")
    include("sfcc_tests.jl")
    include("solver_tests.jl")
    include("inference_tests.jl")
end
