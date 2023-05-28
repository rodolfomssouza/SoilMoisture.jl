using SoilMoisture
using Test

@testset "SoilMoisture.jl" begin
    # Write your tests here.
    @test leakage(0.5, 0.75, 3.5, 200) == 0.0
    @test leakage(1.0, 0.75, 3.5, 200) == 200
end
