using PspGen
using DifferentialEquations
using Plots
using Roots
using IntervalArithmetic, IntervalRootFinding

# # schordinger integration by ode
# r0 = 1e-6
# Z = 1
# l = 0
# ε = -13.6
# xs = r0:0.01:20.0
# ys = - Z ./ xs
# pot = PspGen.Potential(xs, ys)
# R0, R0p = PspGen.bc_origin(r0, Z, l, ε, pot)
# u0 = [R0, R0p]
# rm = 6.0    # 
# rspan = (r0, rm)
# p = [l, ε, pot]

function SE_ODE!()
    # SCHRODINGER 2nd order ODE
    #
    #                    d R
    # y_1 = R     y_2 =  ---
    #                    d r
    #
    # f_1 = y_2
    #
    #         2        l(l+1)
    # f_2 = - - y_2 + [------ + 2(v_ks - e)] y_1
    #         r         r^2       
    # l, ε, pot = p
    # rinv = 1 / r
    # r2i = rinv * rinv
    # pot -> Potential
end

function coulomb_ODE!(du, u, p, r)
    v(r) = -1/r
    f(r, e) = v(r) - e
    du[1] = u[2]
    du[2] = f(r, p) * u[1]
end

"""
The log der diff at turning point rm at given eigenvalue e.
"""
function coulomb_ldd(ε::Float64)
    rm = -1/ε   # f(r, e) = -1/r - e

    prob(u0, rspan, e) = ODEProblem(coulomb_ODE!, u0, rspan, e)

    ODE_INF = 1.0e-12   # The practical value of function at INF.
    
    rmin = 1.0e-8  # the smallest grid x value of potential
    a1 = 1 
    a2 = -1/2
    a3 = 1/12 - ε/6

    r = rmin
    a1 = 1
    a2 = -2
    a3 = 1/3
    a4 = -1/36
    bc = a1 * r + a2 * r^2 + a3 * r^3
    bcp = a1 + 2 * a2 * r + 3 * a3 * r^2 + 4 * a4 * r^3
    bc_zero = [bc, bcp]

    # emax = -1.0e-12 # either the smallest possible eigenvalue or a very small value
    rmax = -log(ODE_INF) / sqrt(-1*ε)   # The solution when set V as 0. use the smallest ε if there are more.

    bc_inf = [ODE_INF, -√(-1/ε) * ODE_INF]

    # outward from origin, and inward from inf
    ldd = PspGen.logder_diff(prob, ε, bc_zero, bc_inf, rmin, rmax, rm)
    ldd
end


ε = find_zero(e->coulomb_ldd(e), (-4.0, -0.21))
println(ε)
