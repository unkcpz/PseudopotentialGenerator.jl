using Test
using PspGen
using QuadGK: quadgk

@testset "HydrogenLike" begin
    # test R_nl
    for n in 1:3, l in 0:n-1, Z in 1:5
        t, residue = quadgk(r -> r^2 * R_nl(n, l, r; Z=4)^2, 0, Inf)
        @test t ≈ 1 atol=1e-8
    end
end