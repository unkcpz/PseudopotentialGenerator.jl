# %%
using PseudopotentialGenerator
using Plots
using ProfileView
gr()

# %%
# Hydrogen-like atom
# Z=2 and 4p
Z = 2
n = 4
l = 1

E_exact = -Z^2 / (2.0 * n^2)
mesh = Mesh(1e-8, 40.0, 1e+7, 8000)
Vf(r) = -Z / r
V = Vf.(mesh.r);
E_window = [-5000.0, -0.0]
E_init = -1000.0

# plot iter 1, 10, 12, 14, 50 (converged)
# %%
# converged line
# E, P, Q = solve_radial_eigenproblem(
@benchmark solve_radial_eigenproblem(
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
    perturb = true,
)

# %%

plot(mesh.r, P, label = "converged", line = (2, :solid), lc = :black)

# %%
using Printf

ctp = 7800
mesh.r[ctp]

Pd = Dict()
Qd = Dict()

E = -0.124

# %%

@benchmark sch_outward(l, Z, E, V, mesh.r, mesh.rp)

# %%
for E in [-0.124, -0.1245, -0.125]
    P_out, Q_out, _ = sch_outward(l, Z, E, V, mesh.r, mesh.rp)
    P_in, Q_in, _ = sch_inward(l, E, V, mesh.r, mesh.rp)

    # align P at ctp
    factor = P_out[ctp] / P_in[ctp]
    P_in *= factor
    Q_in *= factor

    P = zeros(length(mesh.r))
    Q = zeros(length(mesh.r))
    P[1:ctp], P[ctp:end] = P_out[1:ctp], P_in[ctp:end]
    Q[1:ctp], Q[ctp:end] = Q_out[1:ctp], Q_in[ctp:end]

    # normalize
    S = integrate(P .^ 2, mesh.rp, method = :trapz7)
    S = sqrt(abs(S))
    P /= S
    Q /= S

    Pd[@sprintf("%.4f", E)] = P
    Qd[@sprintf("%.4f", E)] = Q
end

p_P = plot(mesh.r, Pd["-0.1240"], label = "E=-0.124(Ha)", line = (2, :dash), lc = :red)
plot!(mesh.r, Pd["-0.1245"], label = "E=-0.1245(Ha)", line = (2, :dashdot), lc = :blue)
plot!(mesh.r, Pd["-0.1250"], label = "E=-0.125(Ha)", line = (2, :solid), lc = :green)
plot!(ylabel = "P(r)")
plot!(title = "Continuity of Z=2, 4p")

p_Q = plot(mesh.r, Qd["-0.1240"], label = "E=-0.124(Ha)", line = (2, :dash), lc = :red)
plot!(mesh.r, Qd["-0.1245"], label = "E=-0.1245(Ha)", line = (2, :dashdot), lc = :blue)
plot!(mesh.r, Qd["-0.1250"], label = "E=-0.125(Ha)", line = (2, :solid), lc = :green)
plot!(xlabel = "r (a.u.)", ylabel = "dP(r)/P(r)")

display(plot(p_P, p_Q, layout = (2, 1)))

# %%
savefig("/tmp/conti_pp.pdf")
