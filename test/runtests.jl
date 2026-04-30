using Test, XLZD, JSON3

@testset "XLZD" begin
    include("test_geometry.jl")
    include("test_physics.jl")
    include("test_mc.jl")
    include("test_output.jl")
end
