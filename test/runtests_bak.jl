using Test
using PspGen
using DifferentialEquations
using Plots

# Atom Be Z=4, n=2
@test PspGen.compute_hydrogen(4, 2) ≈ -2.0 

# Potential get test
p = PspGen.Potential([0.0, 1, 2, 4], [0, 1.1, 2.4, 6])
pflat = PspGen.Potential([0.0, 0.11, 1.0], [1.1, 1.1, 1.1])
@test PspGen.v(1.2, p) ≈ 1.3392
@test PspGen.v(0.5, pflat) ≈ 1.1

# schordinger integration by ode
# r0 = 1e-6
# Z = 1
# l = 2
# ε = -0.4
# pot = PspGen.Potential([0.0, 10.0], [0.0, 0.0])
# R0, R0p = PspGen.bc_origin(r0, Z, l, ε, pot)
# u0 = [R0, R0p]
# rm = 2.0    # 
# rspan = (r0, rm)
# p = [l, ε, pot]
# prob = ODEProblem(PspGen.SE_ODE!, u0, rspan, p)
# sol = solve(prob)
# pl = plot(sol)
# savefig(pl, "myplot.png")

prob(e) = ODEProblem(SE_ODE!, u0, rspan, e)
sol_f(e) = solve(prob(e), Tsit5(), reltol=1e-8, abstol=1e-8)
e_eigen = find_zero(e->sol_f(e).u[end][1], (-1.0, -0.1))
println(e_eigen)
pl = plot(sol_f(e_eigen))
savefig(pl, "dummy_ode.png")


