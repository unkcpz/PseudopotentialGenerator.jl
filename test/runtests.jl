using Test
using PGEN

@testset "ode1d" begin
    x = range(0, 3, length=20)
    y = sin.(x)
    xp = cos.(x)
    for method in [:trapz1, :trapz3, :trapz5, :trapz7, :simpson, :adams]
        f_method = Symbol("f_", method)
        @test integrate(y, xp, method=method) ≈ FPGEN.integrate(y, xp, method=f_method) atol = 1e-12
    end
end

@testset "mesh" begin
    r_min = 1.0
    r_max = 10.0
    a = 1.0
    N = 1
    mesh = FPGEN.mesh(r_min, r_max, a, N)
    println(mesh)
end