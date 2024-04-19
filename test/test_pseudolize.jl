using JLD2

# AE results for different pseudolize methods
# Only need to run once

Z = 6
mesh = Mesh(1e-7, 50.0, 2.7e+6, 10000)
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
    using Plots

    plot(mesh.r, ae_info.vae, label="AE")
    ylims!(-10, 10)
    xlims!(0, 20)

    v_pspot = pseudolize(ae_info, mesh, rc; method=:TM, kbform=false)    

    for nl in keys(rc)
        vline!([rc[nl]], label="rc: l=$(nl.l)", linestyle=:dash)
        plot!(mesh.r, v_pspot.v[nl], label="pspot l=$(nl.l)")
    end

    # save it to file
    savefig("pseudolize_TM.png")

    # TODO: consistency check the ps_pot will give the expected ps_wave
end