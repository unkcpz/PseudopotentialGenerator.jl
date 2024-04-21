@testset "reigen norel" begin
    # Hydrogen-like atom
    Z = 3
    using Plots
    using Colors

    palette = distinguishable_colors(10)  # Generate 10 distinct colors

    plot()
    for perturb in [false], n in 1:2, l in 0:n-1
    #for perturb in [true], n in 1:1, l in 0:n-1
        E_exact = -Z^2 / (2.0 * n^2)
        mesh = Mesh(1e-8, 50.0, 1e+7, 20000)
        V(r) = -Z / r
        E_window = [-5000.0, -0.0]
        E_init = -1000.0

        E, P, Q = solve_radial_eigenproblem(n, l, Z, V, mesh; tol=1e-9, max_iter=100, E_window=E_window, E_ini = E_init, rel = false, perturb = perturb)
        @test E ≈ E_exact atol = 1e-4

        Vinner = @. V(mesh.r) + l*(l+1) / (2 * mesh.r^2)
        ctp = find_ctp(Vinner, E)

        # random color
        color = palette[rand(1:10)]

        plot!(mesh.r, P, label="n = $n, l = $l", color=color)
        # plot vertical line at ctp
        vline!([mesh.r[ctp]], label="ctp: n=$n, l=$l", linestyle=:dash, color=color)
        savefig("test.png")

        #E_fpgen, P_fpgen, Q_fpgen = FPGEN.solve_radial_eigenproblem(n, l, Z, V, mesh; tol=1e-9, max_iter=100, E_window=E_window, E_ini = E_init, rel = false, perturb = perturb)
        #@test E_fpgen ≈ E_exact atol = 1e-4

        # The wavefunction returned should also be the same
        #for idx in [10, 50, 100]
        #    @test Q[idx] ≈ Q_fpgen[idx] atol = 1e-4
        #    @test P[idx] ≈ P_fpgen[idx] atol = 1e-4
        #end
    end
end