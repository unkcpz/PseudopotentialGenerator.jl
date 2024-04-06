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
cached_filename = "pseudo_test_Z_6_orbs.jld2"
if !isfile(cached_filename)
    res = self_consistent_field(Z, mesh, orbs, mixing_beta=0.2, abstol=1e-7, maxiters_scf=200, perturb=true)

    Ys = zeros(Float64, length(orbs), length(mesh.r))
    for (idx, orb) in enumerate(orbs)
        ϕ = res.ϕs[(n=orb.n, l=orb.l)]
        Ys[idx, :] = ϕ
    end
    @save cached_filename res
else
    @load cached_filename res
end

# TODO: rc for every l channels, it can be multiple rc
# will only pseudolize on p orbital for test.
rc = 2.0

# Check integration of ρ to inf is equal to Z
# Check the wavefunction is normalized (∑ϕ^2 dr = Z)
# plot the wavefunction and density

@testset "Pseudolize TM" begin
    using Plots

    vae = res.v_tot
    plot(mesh.r, vae, label="AE")

    # add a vertical dash line to show the rc
    vline!([rc], label="rc", linestyle=:dash)

    for l in 0:1
        nl = (; n=2, l=l)
        εl = res.ε_lst[nl]
        ϕ = res.ϕs[nl]

        v_pspot, ϕ_ps = pseudolize(ϕ, nl, rc, εl, vae, mesh; method=:TM)    

        plot!(mesh.r, v_pspot, label="pspot l=$l")
    end
    ylims!(-10, 10)
    xlims!(0, 20)

    # save it to file
    savefig("pseudolize_TM.png")
end