@testset "mesh" begin
    r_min = 1.0
    r_max = 50.0
    a = 1e+9
    N = 10
    mesh = Mesh(r_min, r_max, a, N)
    @test mesh.r ≈ FPGEN.mesh_exp(r_min, r_max, a, N) atol = 1e-12
    @test mesh.rp ≈ FPGEN.mesh_exp_deriv(r_min, r_max, a, N) atol = 1e-12
    @test mesh_exp_deriv2(mesh) ≈ FPGEN.mesh_exp_deriv2(r_min, r_max, a, N) atol = 1e-12

    a = 1.0 # uniform grid
    mesh = Mesh(r_min, r_max, a, N)
    @test mesh.r ≈ FPGEN.mesh_exp(r_min, r_max, a, N) atol = 1e-12
end