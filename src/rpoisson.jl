using DifferentialEquations

function poisson_outward(ρ::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64})::Vector{Float64}
    # Converted 2nd Poisson ODE in rgrid to uniform grid
    # u1p = u2 * Rp
    # u2p = -(4*pi*rho + 2*u2/r) * Rp

    r_mid = midpoints(r)
    rp_mid = midpoints(rp, r)
    ρ_mid = midpoints(ρ, r)

    function rp_u(idx)
        # idx is integer get from vec otherwise get from midpoints
        if isinteger(idx)
            idx = Int(idx)
            return rp[idx]
        else
            return rp_mid[floor(Int, idx)]
        end
    end

    function r_u(idx)
        # idx is integer get from vec otherwise get from midpoints
        if isinteger(idx)
            idx = Int(idx)
            return r[idx]
        else
            return r_mid[floor(Int, idx)]
        end
    end

    function ρ_u(idx)
        # idx is integer get from vec otherwise get from midpoints
        if isinteger(idx)
            idx = Int(idx)
            return ρ[idx]
        else
            return ρ_mid[floor(Int, idx)]
        end
    end

    function poisson!(du, u, _p, t)
        du[1] = u[2] * rp_u(t)
        du[2] = -(4π * ρ_u(t) + 2 * u[1] / r_u(t)) * rp_u(t)
        nothing
    end

    # step t is integer and t = 1, 2, ..., N
    N = length(r)
    y0 = [4π * integrate(ρ .* r, rp), 0.0]
    prob = DiscreteProblem(poisson!, y0, (1, N))
    sol = solve(prob, ABM54(), dt=1)

    sol[1, :]
end

function rpoisson_outward_pc(ρ::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64})::Vector{Float64}
    # Translate from dftatom:
    # It rewrites it to the equivalent system of first order ODEs on a uniform
    # grid:
    #   u1 = V
    #   u2 = V'
    #   u1p = u2 * Rp
    #   u2p = -(4*pi*rho + 2*u2/r) * Rp
    # and integrates outward using Adams method. The initial conditions are:
    #   V (R(1)) = u1(1) = 4*pi * \int r * rho(r) dr
    #   V'(R(1)) = u2(1) = 0
    N = length(r)

    u1 = zeros(Float64, N)
    u2 = zeros(Float64, N)
    u1p = zeros(Float64, N)
    u2p = zeros(Float64, N)

    if N < 5
        throw(ArgumentError("Too less points to integrate, N must be greater than 5"))
    end

    # Integrate using the Runge-Kutta 4th order method for the first 4 points
    # Boundary condition
    y0 = [4π * integrate(ρ .* r, rp), 0.0]

    C1 = -4π .* ρ[1:4]
    C2 = -2 ./ r[1:4]
    ρ_mid = midpoints(ρ[1:4], r[1:4])
    rmid = midpoints(r[1:4])

    C1mid = -4π .* ρ_mid[1:3]
    C2mid = -2 ./ rmid[1:3]
    y1, y2, _ = rk4_integrate_poisson(r[1:4], y0, C1, C2, C1mid, C2mid, 1e10)
    u1[1:4] .= y1[1:4]
    u2[1:4] .= y2[1:4]
    @. u1p[1:4] = u2[1:4] * rp[1:4]
    @. u2p[1:4] = -(4π * ρ[1:4] + 2 * u2[1:4] / r[1:4]) * rp[1:4]

    # Adams method for the rest of the points
    for i in 4:N-1
        u1[i+1] = u1[i] + adams_extrapolation_outward(u1p, i)
        u2[i+1] = u2[i] + adams_extrapolation_outward(u2p, i)
        for _ in 1:2
            u1p[i+1] = u2[i+1] * rp[i+1]
            u2p[i+1] = -(4π * ρ[i+1] + 2 * u2[i+1] / r[i+1]) * rp[i+1]
            u1[i+1] = u1[i] + adams_interp_outward(u1p, i)
            u2[i+1] = u2[i] + adams_interp_outward(u2p, i)
        end
    end

    u1
end

function rpoisson_outward_pc(ρ::Vector{Float64}, mesh::Mesh)::Vector{Float64}
    u1 = rpoisson_outward_pc(ρ, mesh.r, mesh.rp)

    u1
end