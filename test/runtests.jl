using Test, NJLModels
using NJLModels.Zablocki

param = Zablocki.Parameters2()
T = 0.5 * rand()
mu = 0.5 * rand()

@testset "Gap Equation" begin
    @test [0.0, -0.1] < massgap(T, mu, param).zero < [param.Λ, param.Λ]
    @test 0.0 < Zablocki.massgap_m(T, mu, param) < param.Λ
    @test 0.2 < massgap(0.01, 0.0, param).zero[1] < 0.45
end
