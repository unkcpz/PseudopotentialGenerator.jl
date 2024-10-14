# NOTE: To run this example the FPGEN part require compile wrapper in deps
# Goto deps folder and run `make all`

# %%
using PseudopotentialGenerator
using Plots
using OrdinaryDiffEq
gr()

# %%
# Hydrogen-like atom
# Z=92 and 6d
Z = 92
n = 6
l = 2

E_exact = -Z^2 / (2.0 * n^2)
mesh = Mesh(1e-10, 100.0, 1e+7, 10000)
Vf(r) = -Z / r
V = Vf.(mesh.r);
E_window = [-5000.0, -0.0]
E_init = -1000.0

# plot iter 1, 10, 12, 14, 50 (converged)
j % %
# converged line
E, P_ref, Q_sol = solve_radial_eigenproblem(
    n,
    l,
    Z,
    V,
    mesh;
    tol = 1e-9,
    max_iter = 100,
    E_window = E_window,
    E_ini = E_init,
    rel = false,
    perturb = false,
)

plot(mesh.r, P_ref, label = "out/in-ward Ref", line = (4, :solid), lc = :blue)
plot!(xlim = (0, 4))
plot!(ylim = (-2, 2))

ctp = 5500
mesh.r[ctp]
P_ref[ctp]

E = E_exact

# ABM5
P_out = zeros(length(mesh.r))
Q_out = zeros(length(mesh.r))
P_out[1:end], Q_out[1:end], _ =
    sch_outward(l, Z, E, V[1:end], mesh.r[1:end], mesh.rp[1:end]; algo = VCABM5())

# align P to P_ref at ctp
factor = P_ref[ctp] / P_out[ctp]
P_out *= factor
Q_out *= factor

plot!(xlim = (0, 4))
plot!(ylim = (-2, 2))
plot!(xlabel = "r (a.u.)")
plot!(mesh.r[1:end], P_out[1:end], label = "VCABM5", line = (2, :dashdot))

# ABM4
P_out = zeros(length(mesh.r))
Q_out = zeros(length(mesh.r))
P_out[1:end], Q_out[1:end], _ =
    sch_outward(l, Z, E, V[1:end], mesh.r[1:end], mesh.rp[1:end]; algo = VCABM4())

# align P to P_ref at ctp
factor = P_ref[ctp] / P_out[ctp]
P_out *= factor
Q_out *= factor

plot!(xlim = (0, 4))
plot!(ylim = (-2, 2))
plot!(mesh.r[1:end], P_out[1:end], label = "VCABM4", line = (2, :dashdot))

# ABM54 without adaptive
P_out = zeros(length(mesh.r))
Q_out = zeros(length(mesh.r))
P_out[1:end], Q_out[1:end], _ =
    sch_outward(l, Z, E, V[1:end], mesh.r[1:end], mesh.rp[1:end]; algo = ABM54())

# align P to P_ref at ctp
factor = P_ref[ctp] / P_out[ctp]
P_out *= factor
Q_out *= factor

plot!(xlim = (0, 4))
plot!(ylim = (-2, 2))
plot!(mesh.r[1:end], P_out[1:end], label = "ABM54", line = (2, :dashdot))


P_out = zeros(length(mesh.r))
Q_out = zeros(length(mesh.r))
P_out[1:end], Q_out[1:end], _ =
    FPGEN.schroed_outward_adams(l, Z, E, Vf, mesh.r[1:end], mesh.rp[1:end])

# align P to P_ref at ctp
factor = P_ref[ctp] / P_out[ctp]
P_out *= factor
Q_out *= factor

plot!(mesh.r[1:end], P_out[1:end], label = "DFTATOM (4th order adams)", line = (2, :dash))

# %%
savefig("/tmp/stab_pp.pdf")
