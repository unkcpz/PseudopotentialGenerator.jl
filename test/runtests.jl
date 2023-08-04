using Test
using PspGen
using QuadGK: quadgk

# # lschfb test with woods_saxon 
# u0 = -50.0
# a = 0.6
# r0 = 7.0
# u1 = - u0 / a
# Z = 1
# rgrid, dr = compute_rgrid(Z)
# v_ws = woods_saxon(rgrid, u0, u1, r0, a)

# Z = 1
# n = 1
# l = 0
# alpha = 0.05
# v_hulthen = hulthen(rgrid, Z, alpha)
# sol_orb(n, l, rgrid, v_hulthen, Z, false)


@testset "AE scf" begin
    # test scf! AE
    # Hydrogen
    Z = 1
    rgrid, dr = compute_rgrid(Z)
    srel = false
    xc = (:lda_x, :lda_c_pz)
    orbs = [Orbital(1, 0, 1, Z, rgrid)]   # 1s1

    _, _, E = scf!(Z, orbs, rgrid, dr, xc, srel)
    println("Total energy: ", E)
    # orbs is updated
    l = 0
    rc = 1.54
    nbessel = 3
    aewfc = orbs[1].ur
    # print(aewfc)
    rrkj(aewfc, l, rc, rgrid, dr, nbessel)
end

@testset "HydrogenLike" begin
    t, residue = quadgk(r -> r^2 * R_nl(2, 1, r; Z=4)^2, 0, Inf)
    @test t == 1
    @test residue < 1e-8
end