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

@testset "mesh derivative" begin
    mesh = Mesh(1.0, 50.0, 1e+9, 1000)
    y = sin.(mesh.r)
    expected_yp = cos.(mesh.r)

    for ic in [100, 200, 300, 400, 500]
        yp = dfdr(y, mesh, ic)
        ypp = d2fdr2(y, mesh, ic)
        @test yp ≈ expected_yp[ic] atol = 1e-8
        # @test ypp ≈ -y[ic] atol = 1e-2 # TODO: why d2d not accurate?
    end
end