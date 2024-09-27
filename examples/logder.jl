# %%
using PGEN
using Plots
gr()

# Define different color for each l
colors = Dict(0 => :red, 1 => :blue, 2 => :green, 3 => :orange, 4 => :purple)

# %%
# AE results for different pseudolize methods
# Only need to run once

Z = 6
mesh = Mesh(1e-7, 50.0, 2.7e+6, 10000);
orbs = [Orbital(Z, 1, 0, 2), Orbital(Z, 2, 0, 2), Orbital(Z, 2, 1, 2)];

res = self_consistent_field(
    Z,
    mesh,
    orbs,
    mixing_beta = 0.2,
    abstol = 1e-7,
    maxiters_scf = 200,
    perturb = true,
)

ae_info = (ε_lst = res.ε_lst, ϕs = res.ϕs, vae = res.v_tot, occs = res.occs, xc = res.xc);

# %%

# rc for every orbital, pass as a dict
# angular momentum => rc
rc = Dict{NamedTuple{(:n, :l),Tuple{Int64,Int64}},Float64}(
    (n = 2, l = 0) => 1.9,
    (n = 2, l = 1) => 1.8,
)

# Check integration of ρ to inf is equal to Z
# Check the wavefunction is normalized (∑ϕ^2 dr = Z)
# plot the wavefunction and density
v_pspot, ϕ_ps = pseudolize(ae_info, mesh, rc; method = :TM, kbform = false)


# %%
# Plot arctan logder and compare with AE
# start a new plot

rcx = 2.0 # the r cut to compute atan logder
x_ae = Dict()
y_ae = Dict()
x_ps = Dict()
y_ps = Dict()

for nl in keys(rc)
    l = nl.l
    x_ae[nl], y_ae[nl] =
        compute_atanld(l, Z, mesh, ae_info.vae, rcx, window = [-10.0, 10.0], δ = 0.05)

    x_ps[nl], y_ps[nl] =
        compute_atanld(l, Z, mesh, v_pspot.v[nl], rcx, window = [-10.0, 10.0], δ = 0.05)
end

# %%
plot(title = "logrithmic derivative")
plot!(xlabel = "Energy (Ha)", ylabel = "tan-1 psi' / psi")

for nl in keys(rc)
    l = nl.l
    plot!(x_ae[nl], y_ae[nl], label = "AE: l = $l", line = (2, :solid), color = colors[l])

    plot!(x_ps[nl], y_ps[nl], label = "PS: l = $l", line = (2, :dash), color = colors[l])
end

ylims!(-10, 10)
display(plot!())

# %%
savefig("/tmp/atanlogder_TM_SL.pdf")

# %%
# Plot ae wavefunctions
plot_wfc = plot()
plot!(xlims = (0, 8), ylims = (-2, 2))
plot!(title = "Z=$Z, s/p channels")

nl = (n = 2, l = 0)
vline!([rc[nl]], label = "", line = (1, :dot), lc = colors[nl.l])

ic = findfirst(r -> r > rc[nl], mesh.r) - 1

ϕ_ae = ae_info.ϕs[nl];
ae_rc = ϕ_ae[ic];

plot!(
    mesh.r,
    ϕ_ae,
    label = "AE, n=$(nl.n), l=$(nl.l)",
    line = (3, :solid),
    lc = colors[nl.l],
)

phi_ps = ϕ_ps[nl]
ps_rc = phi_ps[ic]
phi_ps = phi_ps * ae_rc / ps_rc
plot!(
    mesh.r,
    phi_ps,
    label = "PS, n=$(nl.n), l=$(nl.l)",
    line = (2, :dash),
    lc = colors[nl.l],
)

nl = (n = 2, l = 1)
vline!([rc[nl]], label = "", line = (1, :dot), lc = colors[nl.l])

ic = findfirst(r -> r > rc[nl], mesh.r) - 1

ϕ_ae = ae_info.ϕs[nl];
ae_rc = ϕ_ae[ic];

plot!(
    mesh.r,
    ϕ_ae,
    label = "AE, n=$(nl.n), l=$(nl.l)",
    line = (3, :solid),
    lc = colors[nl.l],
)

phi_ps = ϕ_ps[nl]
ps_rc = phi_ps[ic]
phi_ps = phi_ps * ae_rc / ps_rc
plot!(
    mesh.r,
    phi_ps,
    label = "PS, n=$(nl.n), l=$(nl.l)",
    line = (2, :dash),
    lc = colors[nl.l],
)
plot!(ylabel = "wavefunction")

# %%
# potential
plot_pot = plot(mesh.r, ae_info.vae, label = "AE pot", line = (3, :solid), lc = :black)
ylims!(-3, 1)
xlims!(0, 8)

for nl in keys(rc)
    vline!([rc[nl]], label = "l=$(nl.l), rc=$(rc[nl])", line = (1, :dot), lc = colors[nl.l])
    plot!(mesh.r, v_pspot.v[nl], label = "", line = (2, :dash), lc = colors[nl.l])
end
plot!(ylabel = "V(r) (Ha)")
plot!(xlabel = "r (a.u.)")

# %%
plot(plot_wfc, plot_pot, layout = (2, 1))

# %%

# save it to file
savefig("/tmp/pseudolize_TM.pdf")

# %%
# Plot arctan logder and compare with AE
# start a new plot

rcx = 2.0 # the r cut to compute atan logder
x_ae = Dict()
y_ae = Dict()
x_ps = Dict()
y_ps = Dict()

for nl in keys(rc)
    l = nl.l
    x_ae[nl], y_ae[nl] =
        compute_atanld(l, Z, mesh, ae_info.vae, rcx, window = [-10.0, 10.0], δ = 0.05)

    x_ps[nl], y_ps[nl] =
        compute_atanld(l, Z, mesh, v_pspot.v[nl], rcx, window = [-10.0, 10.0], δ = 0.05)
end

# %%
rcx = 2.4
x_ae_large_rc = Dict()
y_ae_large_rc = Dict()
x_ps_large_rc = Dict()
y_ps_large_rc = Dict()
for nl in keys(rc)
    l = nl.l
    x_ae_large_rc[nl], y_ae_large_rc[nl] =
        compute_atanld(l, Z, mesh, ae_info.vae, rcx, window = [-10.0, 10.0], δ = 0.05)

    x_ps_large_rc[nl], y_ps_large_rc[nl] =
        compute_atanld(l, Z, mesh, v_pspot.v[nl], rcx, window = [-10.0, 10.0], δ = 0.05)
end

# %%
p_diff_rcx = plot(title = "(compare between rcx = 2.0 and rcx = 2.4)")
plot!(xlabel = "Energy (Ha)", ylabel = "tan-1 psi' / psi")
ylims!(-12, 10)

for nl in keys(rc)
    l = nl.l
    plot!(
        x_ae[nl],
        y_ae[nl],
        label = "AE: l = $l (rcx=2.0)",
        line = (2, :solid),
        color = colors[l],
    )

    plot!(
        x_ae_large_rc[nl],
        y_ae_large_rc[nl],
        label = "AE rc_large: l = $l (rc=2.4)",
        line = (2, :dash),
        color = colors[l],
    )
end

p_rcx_small = plot(title = "rcx = 2.0: compare between AE and PS")
plot!(xlabel = "Energy (Ha)", ylabel = "tan-1 psi' / psi")
ylims!(-12, 10)

for nl in keys(rc)
    l = nl.l
    plot!(
        x_ae[nl],
        y_ae[nl],
        label = "AE: l = $l (rcx=2.0)",
        line = (2, :solid),
        color = colors[l],
    )

    plot!(
        x_ps[nl],
        y_ps[nl],
        label = "PS: l = $l (rcx=2.0)",
        line = (2, :dash),
        color = colors[l],
    )

    plot!(
        x_ae[nl],
        y_ae[nl] - y_ps[nl] .- 10,
        label = "AE-PS: l = $l (rcx=2.0)",
        line = (2, :dot),
        color = colors[l],
    )
end


p_rcx_large = plot(title = "rcx = 2.4: compare between AE and PS")
plot!(xlabel = "Energy (Ha)", ylabel = "tan-1 psi' / psi")
ylims!(-12, 10)

for nl in keys(rc)
    l = nl.l
    plot!(
        x_ae_large_rc[nl],
        y_ae_large_rc[nl],
        label = "AE: l = $l (rcx=2.4)",
        line = (2, :solid),
        color = colors[l],
    )

    plot!(
        x_ps_large_rc[nl],
        y_ps_large_rc[nl],
        label = "PS: l = $l (rcx=2.4)",
        line = (2, :dash),
        color = colors[l],
    )
    plot!(
        x_ae_large_rc[nl],
        y_ae_large_rc[nl] - y_ps_large_rc[nl] .- 10,
        label = "AE-PS: l = $l (rcx=2.4)",
        line = (2, :dot),
        color = colors[l],
    )
end


display(plot(p_diff_rcx, p_rcx_small, p_rcx_large, layout = (3, 1), size = (600, 900)))

# %%
