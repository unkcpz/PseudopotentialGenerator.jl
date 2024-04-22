using OrdinaryDiffEq

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
    rmids = midpoints(r)
    Vmids = midpoints(V, r)
    Cmids = @. 2 * (Vmids - E) + l*(l+1) / rmids^2
    rpmids = r[2:end] - r[1:end-1]
    # TRY??: rp_mid = sqrt.(rp[1:end-1] .* rp[2:end])    # Only for exponential mesh, compare with rp_mid above, err < 1e-10, not help

    # u1p = u2 * rp
    # u2p = C * u1 * rp
    function f!(du, u, _p, t)
        _t = floor(Int, t)
        if _t == t
            rpt = rp[_t]
            Ct = C[_t]
        else
            rpt = rpmids[_t]
            Ct = Cmids[_t]
        end
        du[1] = u[2] * rpt
        du[2] = Ct * u[1] * rpt
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

    rmids = midpoints(r)
    Vmids = midpoints(V, r)
    Cmids = @. 2 * (Vmids - E) + l*(l+1) / rmids^2
    rpmids = r[2:end] - r[1:end-1]
    # TRY??: rp_mid = sqrt.(rp[1:end-1] .* rp[2:end])    # Only for exponential mesh, compare with rp_mid above, err < 1e-10, not help

    # u1p = u2 * rp
    # u2p = C * u1 * rp
    function f!(du, u, _p, t)
        _t = floor(Int, t)
        if _t == t
            rpt = rp[_t]
            Ct = C[_t]
        else
            rpt = rpmids[_t]
            Ct = Cmids[_t]
        end
        du[1] = u[2] * rpt
        du[2] = Ct * u[1] * rpt
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
