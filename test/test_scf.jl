@testset "NL solver mixing" begin
    # Test the mixing NL solver
    Z = 5
    mesh = Mesh(1e-7, 50.0, 2.7e+6, 6000)

    function regular_atom_orbitals(Z)::Vector{Orbital}
        if Z == 1
            orbs = [
                Orbital(Z, 1, 0, 1),
            ]
        end

        if Z == 3
            orbs = [
                Orbital(Z, 1, 0, 2),
                Orbital(Z, 2, 0, 1),
            ]
        end

        if Z == 5
            orbs = [
                Orbital(Z, 1, 0, 2),
                Orbital(Z, 2, 0, 2),
                Orbital(Z, 2, 1, 1),
            ]
        end

        if Z == 6
            orbs = [
                Orbital(Z, 1, 0, 2),
                Orbital(Z, 2, 0, 2),
                Orbital(Z, 2, 1, 2),
            ]
        end
        orbs
    end

    orbs = regular_atom_orbitals(Z)
    res = self_consistent_field(Z, mesh, orbs, mixing_beta=0.2, abstol=1e-7, maxiters_scf=200, perturb=true)

    # Reference data from:
    # https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations
    @test res.E_tot ≈ -24.344198 atol = 1e-6

    #using Plots
    #ϕ_1s = res.ϕs[(n=1, l=0)]
    #ϕ_2s = res.ϕs[(n=2, l=0)]
    #ϕ_2p = res.ϕs[(n=2, l=1)]

    #plot(mesh.r, ϕ_1s, label="1s")
    #plot!(mesh.r, ϕ_2s, label="2s")
    #plot!(mesh.r, ϕ_2p, label="2p")

    #savefig("nl_solver_mixing.png")
end