using Test
using PGEN

#include("test_mesh.jl")
#include("test_scf.jl")
#include("test_integ.jl")
include("test_eigensolver.jl")
#include("test_pseudolize.jl")

#using Plots
#
#mesh = Mesh(1e-6, 50.0, 1e+9, 5000)
#
#Z = 1
#n = 1
#l = 0
#E_exact = -Z^2 / (2.0 * n^2)
#V = -Z ./ mesh.r
#E_window = [-5000.0, -0.0]
#E_init = -1000.0
#E = E_exact
#N = length(mesh.r)
#
##E, P, Q = solve_radial_eigenproblem(n, l, Z, V, mesh; tol=1e-9, max_iter=100, E_window=E_window, E_ini = E_init, rel = false, perturb = true, normalize=false)
#
#P, Q, imax = schroed_inward_adams(l, E, V, mesh.r, mesh.rp)
#P_, Q_, _ = sch_inward(l, E, V, mesh.r, mesh.rp)
#
#println("E = ", E)
#
#plot()
#plot(mesh.r, P, label="P")
#plot!(mesh.r, P_, label="P_")
#plot!(mesh.r, P-P_, label="P-P_")
#
#plot!(mesh.r, Q, label="Q")
#plot!(mesh.r, Q_, label="Q_")
#plot!(mesh.r, Q-Q_, label="Q-Q_")
#xlims!(0, 50)
#ylims!(-0.5, 2)
#
#savefig("test.png")
#
##Q, P, imin = schroed_inward_adams(l, E, V, mesh.r, mesh.rp)
#
#
#
#