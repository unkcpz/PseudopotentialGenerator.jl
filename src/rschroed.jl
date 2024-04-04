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
    u1 = zeros(Float64, N)  # Q
    u2 = zeros(Float64, N)  # P
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
    u1 = zeros(Float64, N)  # Q
    u2 = zeros(Float64, N)  # P
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
        # throw(ArgumentError("No enough starting points to integrate inward"))
        imax = N-5
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