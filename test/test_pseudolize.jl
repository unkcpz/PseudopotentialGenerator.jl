using JLD2
using Plots

# AE results for different pseudolize methods
# Only need to run once

Z = 6
mesh = Mesh(1e-7, 20.0, 2.7e+6, 1000)
orbs = [
    Orbital(Z, 1, 0, 2),
    Orbital(Z, 2, 0, 2),
    Orbital(Z, 2, 1, 2),
]

# caching the result of array of each orbital to jls files
# named with pseudo_test_Z_6_orb_<n>_<l>.jld2
cached_filename = "pseudo_test_Z_6_ae_info.jld2"
if !isfile(cached_filename)
    res = self_consistent_field(Z, mesh, orbs, mixing_beta=0.2, abstol=1e-7, maxiters_scf=200, perturb=true)

    ae_info = (
        ε_lst = res.ε_lst,
        ϕs = res.ϕs,
        vae = res.v_tot,
        occs = res.occs,
        xc = res.xc,
    )
    @save cached_filename ae_info
else
    @load cached_filename ae_info
end

# rc for every orbital, pass as a dict
# angular momentum => rc
rc = Dict{NamedTuple{(:n, :l), Tuple{Int64, Int64}}, Float64}(
    (n=2, l=0) => 2.0,
    (n=2, l=1) => 1.9
)

# Check integration of ρ to inf is equal to Z
# Check the wavefunction is normalized (∑ϕ^2 dr = Z)
# plot the wavefunction and density

@testset "Pseudolize TM" begin

    plot(mesh.r, ae_info.vae, label="AE")
    ylims!(-10, 10)
    xlims!(0, 20)

    v_pspot, _ = pseudolize(ae_info, mesh, rc; method=:TM, kbform=false)    

    for nl in keys(rc)
        vline!([rc[nl]], label="rc: l=$(nl.l)", linestyle=:dash)
        plot!(mesh.r, v_pspot.v[nl], label="pspot l=$(nl.l)")
    end

    # save it to file
    savefig("pseudolize_TM.png")

    # Plot arctan logder and compare with AE
    # start a new plot
    plot()

    # Define different color for each l
    colors = Dict(
        0 => :red,
        1 => :blue,
        2 => :green,
        3 => :orange,
        4 => :purple,
    )

    rcx = 2.0 # the r cut to compute atan logder
    for nl in keys(rc)
        l = nl.l
        x_ae, y_ae = compute_atanld(l, Z, mesh, ae_info.vae, rcx, window=[-5.0, 5.0], δ=0.1)
        plot!(x_ae, y_ae, label="AE: l = $l", linestyle=:dash, color=colors[l])

        x_ps, y_ps = compute_atanld(l, Z, mesh, v_pspot.v[nl], rcx, window=[-5.0, 5.0], δ=0.1)
        plot!(x_ps, y_ps, label="PS: l = $l", color=colors[l])
    end

    ylims!(-10, 10)
    savefig("atanlogder_TM_SL.png")
end
