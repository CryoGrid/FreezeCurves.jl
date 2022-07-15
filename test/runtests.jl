using FreezeCurves
using Test
using Unitful

@testset "FreezeCurves.jl" begin
    include("swrc_tests.jl")
    include("sfcc_tests.jl")
    include("solver_tests.jl")
end
