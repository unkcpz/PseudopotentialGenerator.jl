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
        orbs
    end

    orbs = regular_atom_orbitals(Z)
    info = self_consistent_field(Z, mesh, orbs, mixing_beta=0.2, abstol=1e-7, maxiters_scf=200, perturb=true)

    # Reference data from:
    # https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations
    @test info.E_tot ≈ -24.344198 atol = 1e-6
end