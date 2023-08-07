using Test
using PspGen
using QuadGK: quadgk
using Plots
using Dierckx

@testset "Hartree potential and energy" begin
    # Grid
    R = RGrid(3, dx=8e-3, xmin=-9, rmax=100)

    # construct hydrogen_like 3p Li density Z = 3
    u(r) = r * R_nl(3, 1, r; Z=3)
    n(r) = abs2(u(r))

    #t_N, _ = quadgk(r -> n(r), 0, Inf)
    ρ = [n(r) for r in R.r]
    t_N = primitive(R, ρ)
    @test t_N ≈ 1 atol=1e-8

    # test Hartree potential
    ϕ = PspGen.ONCVPSP.hartree(1.0, R.r, ρ)

    # f_sin = sin.(R.r)
    #fd = Spline1D(R.r, f_sin, k=5, bc="nearest")
    #fd1 = derivative(fd, R.r)

    grad_ϕ = gradient(R, ϕ)
    hess_ϕ = hessian(R, ϕ)
    lhs = hess_ϕ + 2 * grad_ϕ ./ R.r 

    xmin = R.r[1]
    #plot(R.r, cos.(R.r), label="A")
    #plot!(R.r, fd1, label="B")
    #plot(R.r, cos.(R.r) - fd1, label="A")
    plot(R.r, lhs, label="LHS", ylims=(-1, 0.3), xlims=(xmin, 0.1))
    plot!(R.r, ρ, label="ρ", ylims=(-1, 1), xlims=(0, 10))
    plot!(R.r, lhs + ρ, label="numerical unstability", ylims=(-1, 1), xlims=(0, 0.04))
    savefig("hartree.png")
end