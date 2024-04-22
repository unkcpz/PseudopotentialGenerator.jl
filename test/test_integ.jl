@testset "Newton-Cotes" begin
    x = range(0, 3, length=20)
    y = sin.(x)
    xp = cos.(x)
    for method in [:trapz1, :trapz3, :trapz5, :trapz7, :simpson]
        f_method = Symbol("f_", method)
        @test integrate(y, xp, method=method) ≈ FPGEN.integrate(y, xp, method=f_method) atol = 1e-12
    end
end

@testset "ode_poisson" begin
    # Test the Poisson solver
    mesh = Mesh(1e-7, 20.0, 1e+5, 1000)

    # wave function in infinite square well potential
    a = 5
    n = 2
    ψ = sqrt(1/a) * sin.(n * π * mesh.r / a)

    # charge density
    ρ = abs2.(ψ)

    # Poisson solver
    vh_expected = FPGEN.rpoisson_outward_pc(ρ, mesh.r, mesh.rp)
    vh = poisson_outward(ρ, mesh)
    @test vh[1:10] ≈ vh_expected[1:10] atol = 1e-10
end