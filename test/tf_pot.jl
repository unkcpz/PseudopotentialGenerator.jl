using Test
using PspGen

@testset "Thomas-Fermi potential" begin
    for Z in 1:2, r in 0.1:0.1:1.0
        @test PspGen.tf(Z, r) ≈ PspGen.ONCVPSP.tf(Z, r)
    end
end