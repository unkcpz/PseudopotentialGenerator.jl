using Test
using PspGen

@testset "radial_grids.jl" begin
    for Z in 1:5
        R = RGrid(Z, dx=1e-5, xmin=-8, rmax=10)
        @test diff(R.r) ≈ R.dr[1:end-1] atol=1e-4
    end
end