using InterfaceMechanics
using Test

# @test hello("Julia") == "Hello, Julia"
# @test domath(2.0) ≈ 7.0

@testset "Test InterfaceMechanics.jl" begin
    @testset "test Couloumb interface" begin
        include("test_couloumb.jl")
    end
end
