using Test
using PspGen
using Plots

@testset "radial_grids.jl" begin
    for Z in 1:5
        R = RGrid(Z, dx=1e-5, xmin=-8, rmax=10)
        @test diff(R.r) ≈ R.dr[1:end-1] atol=1e-4
    end

    # test primitive on hydrogen_like
    for n in 1:3, l in 0:n-1, Z in 1:5
        R = RGrid(Z, dx=8e-3, xmin=-8, rmax=100)
        u(r) = r * R_nl(n, l, r; Z=Z)
        ρ = [abs2(u(r)) for r in R.r]
        t_N = primitive(R, ρ)
        @test t_N ≈ 1 atol=1e-8
    end

    # test numerical derivative
    R = RGrid(1, xmin=-8, rmax=100)
    f_sin = sin.(R.r)
    hess_f = hessian(R, f_sin)

    # Not stable for large r
    plot(R.r, f_sin + hess_f, label="sin(x) - sin''(x)")
    savefig("radial_opt.png")
end