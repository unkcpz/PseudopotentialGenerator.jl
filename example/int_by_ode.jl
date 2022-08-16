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

# function columbode!(du,u,p,r)
#     v(r)=-1/r
#     ri = 1/r
#     r2i = ri * ri
#     l = 0
#     du[1]=u[2]
#     # du[2]=- 2 * ri * u[2] + (l * (l + 1) * r2i + (v(r) - p))*u[1]
#     du[2]= (v(r) - p)*u[1]
# end

# function plotf(rmax)
#     u0=[0.0,1.0]
#     rspan=(1e-10,rmax)
#     prob(e)=ODEProblem(columbode!,u0,rspan,e)
#     solf(e)=solve(prob(e), Tsit5(), reltol=1e-8, abstol=1e-8)
#     e_eigen=find_zero(e->solf(e).u[end][1],(-1.0,-0.1))
#     println(e_eigen)
    
#     solf(e_eigen)
# end

# sol = plotf(30)
# # pl = plot(sol.t, [u[1] for u in sol.u])
# pl = plot(sol)
# savefig(pl, "myplot.png")

# u0=[0.0,1.0]
# rmax = 50.0
# rspan=(1e-10,rmax)
# prob(e) = ODEProblem(PspGen.SE_ODE!, u0, rspan, e)
# sol_f(e) = solve(prob(e), Tsit5(), reltol=1e-8, abstol=1e-8)
# e_eigen = find_zero(e->sol_f(e).u[end][1], (-1.0, -0.1))
# println(e_eigen)
# pl = plot(sol_f(e_eigen))
# savefig(pl, "dummy_ode.png")

ε = find_zero(e->PspGen.ldd_rm(e), (-4.0, -0.21))
println(ε)

# # plot the result function
# ODE_INF = 1.0e-20
# rmax = -log(ODE_INF) / sqrt(-1*ε)

# rmin = 1.0e-10  # the smallest grid x value of potential
# a1 = 1 
# a2 = -1/2
# a3 = 1/12 - ε/6

# r = rmin
# a1 = 1
# a2 = -2
# a3 = 1/3
# a4 = -1/36
# bc = a1 * r + a2 * r^2 + a3 * r^3
# bcp = a1 + 2 * a2 * r + 3 * a3 * r^2 + 4 * a4 * r^3
# bc_zero = [bc, bcp]
# # println(PspGen.ldd_rm(-0.2500))
# prob(u0, rspan, e) = ODEProblem(PspGen.SE_ODE!, u0, rspan, e)
# pl = plot(solve(prob(bc_zero, [rmin, rmax/2], ε), Tsit5(), reltol=1e-12, abstol=1e-12))
# savefig(pl, "dummy_ode.png")