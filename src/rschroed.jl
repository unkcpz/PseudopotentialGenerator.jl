using DifferentialEquations

function sch_outward(l::Int64, Z::Int64, E::Float64, V::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+6)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    C = @. 2 * (V - E) + l*(l+1) / r^2

    # Boundary condition
    # P0 = r ** (l + 1)
    # Q0 = (l + 1) * r ** l
    rmin = r[1]
    if l == 0
        y0 = 1 - Z * rmin
        yp0 = -Z
    else
        y0 = rmin ^ l
        yp0 = l * rmin ^ (l - 1)
    end
    u0 = [y0 * rmin, yp0 * rmin + y0]

    # Compare to dftatom RK4+adams ode solver, this still not stable at large r when rmax is 
    # too large, for Z=92 the coulomb potential is fine with rmax=20, but for rmax=50, the 
    # solution diverged. The reason maybe the interval is too large while the interpolation 
    # is not accurate enough.
    # TODO: rethink how to precisely define the "mid" points
    r_mid = midpoints(r)
    rp_mid = r[2:end] - r[1:end-1]
    #rp_mid = sqrt.(rp[1:end-1] .* rp[2:end])    # Only for exponential mesh, compare with rp_mid above, err < 1e-10, not help
    V_mid = midpoints(V, r)
    C_mid = @. 2 * (V_mid - E) + l*(l+1) / r_mid^2

    # u1p = u2 * rp
    # u2p = C * u1 * rp
    function f!(du, u, _p, t)
        _t = floor(Int, t)
        if _t == t
            du[1] = u[2] * rp[_t]
            du[2] = C[_t] * u[1] * rp[_t]
        else
            du[1] = u[2] * rp_mid[_t]
            du[2] = C_mid[_t] * u[1] * rp_mid[_t]
        end
        nothing
    end

    N = length(r)

    function condition(u, t, integrator)
        abs(u[1]) > max_val || abs(u[2]) > max_val
    end

    function affect!(integrator)
        terminate!(integrator)
    end

    cb = DiscreteCallback(condition, affect!)

    prob = DiscreteProblem(f!, u0, (1, N))
    sol = solve(prob, ABM54(), dt=1, adaptive=false, callback=cb)

    P = zeros(Float64, N)
    Q = zeros(Float64, N)

    imax = length(sol[1, :])

    P[1:imax] .= sol[1, :]
    Q[1:imax] .= sol[2, :]

    P, Q, imax
end

function sch_inward(l::Int64, E::Float64, V::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+6)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    C = @. 2 * (V - E) + l*(l+1) / r^2

    # For large r, the asymptotic is:
    # P(r) = exp(-χ * r)
    # Q(r) = -χ * P(R)
    # where χ = √(-2 * E)
    χ = sqrt(-2 * E)

    # The rmax should not be too large to avoid numerical instability
    # but should be large enough to reach the asymptotic region
    # We make it by apply: exp(-χ * (rmax-rmin)) ~ eps(Float64)
    rmax = r[1] - log(eps(Float64)) / χ

    # find the start idx for the inward integration
    imax = findfirst(r .> rmax)
    if imax === nothing
        # This means the start points may not be enough
        # Use the last points anyway to start the integration
        imax = length(r) - 1
    else
        imax -= 1
    end

    # Boundary condition
    u1 = exp(-χ * r[imax] + χ * r[1])
    u2 = -χ * u1
    u0 = [u1, u2]

    r_mid = midpoints(r)
    rp_mid = r[2:end] .- r[1:end-1]
    V_mid = midpoints(V, r)
    C_mid = @. 2 * (V_mid - E) + l*(l+1) / r_mid^2

    # u1p = u2 * rp
    # u2p = C * u1 * rp
    function f!(du, u, _p, t)
        _t = floor(Int, t)
        if _t == t
            du[1] = u[2] * rp[_t]
            du[2] = C[_t] * u[1] * rp[_t]
        else
            du[1] = u[2] * rp_mid[_t]
            du[2] = C_mid[_t] * u[1] * rp_mid[_t]
        end
        nothing
    end

    function condition(u, t, integrator)
        abs(u[1]) > max_val || abs(u[2]) > max_val
    end

    function affect!(integrator)
        terminate!(integrator)
    end

    cb = DiscreteCallback(condition, affect!)

    prob = DiscreteProblem(f!, u0, (imax, 1))
    sol = solve(prob, ABM54(), dt=-1, adaptive=false, callback=cb)

    P = zeros(Float64, length(r))
    Q = zeros(Float64, length(r))

    n = length(sol[1, :])
    imin = imax - n + 1

    # the integration is from max -> ctp, so we need to reverse the result
    P[imin:imax] .= reverse(sol[1, :])
    Q[imin:imax] .= reverse(sol[2, :])

    P, Q, imin
end


"""
    schroed_outward_adams(l, Z, E, V, r, rp; max_val=1e+6)

Integrate the Schrödinger equation outward using the Adams method.
The Adams method can be applied since the equation is map to the uniform mesh.
"""
function schroed_outward_adams(l::Int64, Z::Int64, E::Float64, V::Function, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+6)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    V = V.(r)
    u1, u2, imax = schroed_outward_adams(l, Z, E, V, r, rp; max_val=max_val)
    u1, u2, imax
end

"""
    schroed_outward_adams(l, Z, E, V, r, rp; max_val=1e+6)

Integrate the Schrödinger equation outward using the Adams method.
The Adams method can be applied since the equation is map to the uniform mesh.
V is a vector.
"""
function schroed_outward_adams(l::Int64, Z::Int64, E::Float64, V::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+6)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    u1 = zeros(Float64, N)  # P
    u2 = zeros(Float64, N)  # Q
    u1p = zeros(Float64, N)
    u2p = zeros(Float64, N)

    if N < 5
        throw(ArgumentError("Too less points to integrate, N must be greater than 5"))
    end

    C = @. 2 * (V - E) + l*(l+1) / r^2

    # Integrate using the Runge-Kutta 4th order method for the first 4 points

    # Boundary condition
    rmin = r[1]
    if l == 0
        y0 = [1 - Z * rmin, -Z]
    else
        y0 = [rmin ^ l, l * rmin ^ (l - 1)]
    end

    # y1 => R(r); y2 => R'(r)
    # P(r) = R(r) * r; Q(r) = R'(r) * r + R(r)
    C1 = @. C[1:4]
    C2 = @. -2 / r[1:4]
    rmid = midpoints(r[1:4])
    Vmid = midpoints(V[1:4], r[1:4])
    C1mid = @. 2 * (Vmid - E) + l*(l+1) / rmid^2
    C2mid = @. -2 / rmid
    y1, y2, _ = rk4_integrate(r[1:4], y0, C1, C2, C1mid, C2mid, max_val)
    @. u1[1:4] = y1[1:4] * r[1:4]
    @. u2[1:4] = y2[1:4] * r[1:4] + y1[1:4]

    # Integrate using the Adams method for the rest of the points
    @. u1p[1:4] = rp[1:4] * u2[1:4]
    @. u2p[1:4] = rp[1:4] * C[1:4] * u1[1:4]

    imax = N
    for i in 4:N-1
        u1p[i] = rp[i] * u2[i]
        u2p[i] = rp[i] * C[i] * u1[i]
        u1_tmp = u1[i] + adams_interp_outward_implicit(u1p, i)
        u2_tmp = u2[i] + adams_interp_outward_implicit(u2p, i)

        λ = 9.0 / 24.0
        Δ = 1 - λ^2 * C[i+1] * rp[i+1]^2
        M11 = 1 / Δ
        M12 = λ * rp[i+1] / Δ
        M21 = λ * C[i+1] * rp[i+1] / Δ
        M22 = 1 / Δ
        M = [M11 M12; M21 M22]

        (u1[i+1], u2[i+1]) = M * [u1_tmp; u2_tmp]

        if abs(u1[i+1]) > max_val || abs(u2[i+1]) > max_val
            imax = i
            break
        end
    end

    u1, u2, imax
end


function schroed_inward_adams(l::Int64, E::Float64, V::Function, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+300)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    V = V.(r)
    u1, u2, imin = schroed_inward_adams(l, E, V, r, rp; max_val=max_val)
    u1, u2, imin
end
"""
    schroed_inward_adams(l, E, V, r, rp; max_val=1e+300)

Integrate the Schrödinger equation inward using the Adams method.
"""
function schroed_inward_adams(l::Int64, E::Float64, V::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64}; max_val::Float64=1e+300)::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    u1 = zeros(Float64, N)  # P
    u2 = zeros(Float64, N)  # Q
    u1p = zeros(Float64, N)
    u2p = zeros(Float64, N)

    if N < 5
        throw(ArgumentError("Too less points to integrate, N must be greater than 5"))
    end

    C = @. 2 * (V - E) + l*(l+1) / r^2

    # For large r, the asymptotic is:
    # P(r) = exp(-χ * r)
    # Q(r) = -χ * P(R)
    # where χ = √(-2 * E)
    χ = sqrt(-2 * E)

    # The rmax should not be too large to avoid numerical instability
    # but should be large enough to reach the asymptotic region
    # We make it by apply: exp(-χ * (rmax-rmin)) ~ eps(Float64)
    rmax = r[1] - log(eps(Float64)) / χ

    # find the start idx for the inward integration
    imax = findfirst(r .> rmax)
    if imax === nothing || imax - 1 > N - 4
        # This means the start points may not be enough
        # But for now just ignore it and start from the end 4 points
        imax = N-4
    else
        imax -= 1
    end

    @. u1[imax:imax+4] = exp(-χ * r[imax:imax+4] + χ * r[1])
    @. u2[imax:imax+4] = -χ * u1[imax:imax+4]
    @. u1p[imax:imax+4] = rp[imax:imax+4] * u2[imax:imax+4]
    @. u2p[imax:imax+4] = rp[imax:imax+4] * C[imax:imax+4] * u1[imax:imax+4]

    imin = 1
    for i in imax:-1:2
        u1p[i] = rp[i] * u2[i]
        u2p[i] = rp[i] * C[i] * u1[i]
        u1_tmp = u1[i] + adams_interp_inward_implicit(u1p, i)
        u2_tmp = u2[i] + adams_interp_inward_implicit(u2p, i)

        λ = - 9.0 / 24.0
        Δ = 1 - λ^2 * C[i-1] * rp[i-1]^2
        M11 = 1 / Δ
        M12 = λ * rp[i-1] / Δ
        M21 = λ * C[i-1] * rp[i-1] / Δ
        M22 = 1 / Δ
        M = [M11 M12; M21 M22]

        (u1[i-1], u2[i-1]) = M * [u1_tmp; u2_tmp]

        if abs(u1[i-1]) > max_val || abs(u2[i-1]) > max_val
            # inward diverge, stop
            imin = i
            break
        end
    end

    u1, u2, imin
end