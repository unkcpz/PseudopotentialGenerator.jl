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
    #@test t_N ≈ 1 atol=1e-8

    # test Hartree potential
    ϕ = hartree(R, ρ, 1.0)

    # Possion
    grad_ϕ = gradient(R, ϕ)
    hess_ϕ = hessian(R, ϕ)
    lhs = hess_ϕ + 2 * grad_ϕ ./ R.r 

    # println(PspGen.find_idx(R, 2.0))
    a = lhs ./ ρ
    # println(a[1300:1400])
    # 0.05

    # Test the stability by compute ρ using Possion Eq
    plot(R.r, lhs, label="LHS", ylims=(-1, 0.3), xlims=(0, 0.1))
    plot!(R.r, ρ, label="ρ", ylims=(-1, 1), xlims=(0, 10))
    plot!(R.r, lhs + ρ, label="numerical unstability", ylims=(-1, 1), xlims=(0, 10))
    savefig("hartree.png")

    # Plot and compare the potential with ONCV
    ϕ_ov = PspGen.ONCVPSP.hartree(1.0, R.r, ρ)
    plot(R.r, ϕ, label="hartree-PE")
    plot!(R.r, ϕ_ov, label="hartree-ONCV", ylims=(-1, 1), xlims=(0, 100))
    savefig("hartree_pot.png")
end