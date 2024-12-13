using OrdinaryDiffEq

function poisson_outward(ρ::Vector{Float64}, mesh::Mesh; max_val = 1e+6, algo = VCABM5())::Vector{Float64}
    # Converted 2nd Poisson ODE in rgrid to uniform grid
    # u1p = u2 * Rp
    # u2p = -(4*pi*rho + 2*u2/r) * Rp

    r = mesh.r
    rp = mesh.rp

    rmids = midpoints(r[1:6])
    rpmids = r[2:end] - r[1:(end - 1)]
    ρmids = midpoints(ρ[1:6], r[1:6])

    function f!(du, u, _p, t)
        _t = floor(Int, t)
        if _t == t
            rpt = rp[_t]
            rt = r[_t]
            ρt = ρ[_t]
        else
            rpt = rpmids[_t]
            rt = rmids[_t]
            ρt = ρmids[_t]
        end
        du[1] = u[2] * rpt
        du[2] = -(4π * ρt + 2 * u[2] / rt) * rpt
        nothing
    end

    N = length(r)

    # Boundary condition
    y0 = [4π * integrate(ρ .* r, rp), 0.0]

    function condition(u, t, integrator)
        abs(u[1]) > max_val || abs(u[2]) > max_val
    end

    function affect!(integrator)
        terminate!(integrator)
    end

    cb = DiscreteCallback(condition, affect!)
    prob = DiscreteProblem(f!, y0, (1, N))
    sol = solve(prob, algo, dt = 1, adaptive = false, callback = cb)

    sol[1, :]
end
