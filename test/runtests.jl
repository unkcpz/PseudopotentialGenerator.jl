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
    r_max = 50.0
    a = 1e+9
    N = 10
    @test mesh_exp(r_min, r_max, a, N) ≈ FPGEN.mesh_exp(r_min, r_max, a, N) atol = 1e-12
    @test mesh_exp_deriv(r_min, r_max, a, N) ≈ FPGEN.mesh_exp_deriv(r_min, r_max, a, N) atol = 1e-12
    @test mesh_exp_deriv2(r_min, r_max, a, N) ≈ FPGEN.mesh_exp_deriv2(r_min, r_max, a, N) atol = 1e-12
end