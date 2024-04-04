function solve_radial_eigenproblem(n::Int64, l::Int64, Z::Int64, V::Function, mesh::Mesh; 
    tol=1e-9, max_iter=200, E_window=[-8000.0, -0.0], E_ini = -1000.0, rel = false, perturb = false)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Bool}
    V = V.(mesh.r)
    E, Q, P, is_converged = solve_radial_eigenproblem(n, l, Z, V, mesh; tol=tol, max_iter=max_iter, E_window=E_window, E_ini = E_ini, rel = rel, perturb = perturb)

    E, Q, P, is_converged
end
function solve_radial_eigenproblem(n::Int64, l::Int64, Z::Int64, V::Vector{Float64}, mesh::Mesh; 
    tol=1e-9, max_iter=200, E_window=[-8000.0, -0.0], E_ini = -1000.0, rel = false, perturb = false)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Bool}
    # TODO: add dirac solver
    N = length(mesh.r)
    P = zeros(Float64, N)
    Q = zeros(Float64, N)
    c = SPEED_OF_LIGHT
    is_converged = false

    E = E_ini

    if n < 1 || l < 0 || l >= n
        throw(ArgumentError("Invalid quantum numbers n=$n, l=$l, l must be in [0, n-1]"))
    end

    Emin, Emax = E_window
    if (E > Emax) || (E < Emin)
        E = (Emin + Emax) / 2
    end

    # shooting method
    last_bisect = true
    for iter in 1:max_iter
        # TODO: write iteration log information
        # bisection converged
        if abs(Emax - Emin) < tol
            if !last_bisect
                # The perturbation theory correction was used in the last
                # iteration and in that case, the consistent stopping criterion is
                # to converge with abs(dE), not abs(Emax - Emin).
                # As such we fail, because abs(Emax - Emin) is converged, but
                # abs(dE) isn't.
                @warn "Bisection converged, but perturbation theory correction was used in the last iteration. The consistent stopping criterion is to converge with abs(dE), not abs(Emax - Emin)."
                break
            end
            
            if abs(Emin - E_window[1]) < eps(Float64)
                # Emin did not change in the loop, something is wrong
                @warn "Emin did not change in the loop, something is wrong"
                break
            end

            if abs(Emax - E_window[2]) < eps(Float64)
                # Emax did not change in the loop, something is wrong
                @warn "Emax did not change in the loop, something is wrong"
                break
            end

            # check if the wave function is not monotonic
            # if it is monotonic, it has no peak
            P, Q, imax = schroed_outward_adams(l, Z, E, V, mesh.r, mesh.rp)
            minidx = get_min_idx(P[1:imax])
            #if is_monotonic(P[1:imax])
            if  minidx <=0
                @warn "Wave function get from outward integration is monotonic, it has no peak"
                break
            end
            P[minidx:end] .= 0.0
            Q[minidx:end] .= 0.0
            
            if count_nodes(P[1:minidx-1]) != n - l - 1
                @warn "Number of nodes in the wave function is not equal to n - l - 1"
                break
            end

            # converged
            # TODO: post-process the wave function
            break
        end

        # Beyond the turning point, meaning the partical is classically forbidden
        # The wave function exponentially decays. 
        Vp = @. V + l*(l+1) / (2 * mesh.r^2)
        ctp = find_ctp(Vp, E)

        if !perturb
            ctp = N
        elseif mesh.r[N] < 2 * mesh.r[ctp] || N < ctp - 10
            # The turning point is too close to the boundary, we need to increase the mesh
            ctp = N
        elseif E > 0
            # The energy is positive, the wave function is oscillatory
            # cannot use inward integration
            ctp = N
        end

        # outward integration
        P[1:ctp], Q[1:ctp], imax = schroed_outward_adams(l, Z, E, V[1:ctp], mesh.r[1:ctp], mesh.rp[1:ctp])
        nnodes = count_nodes(P[1:imax])

        if nnodes != n - l - 1 || ctp == N || imax < ctp
            # FROM DFTATOM: 
            # If the number of nodes is not correct, or we didn't manage to
            # integrate all the way to "ctp", or if "ctp" was too large, we just
            # use bisection:

            # If too many nodes means that the energy is too high
            # otherwise, the energy is too low
            if nnodes > n - l - 1
                Emax = E
            else
                Emin = E
            end
            E = (Emax + Emin) / 2
                     
            last_bisect = true
            continue
        end

        # TODO: The following part should be able to bundle as pertubation method
        # inward integration
        P_inward = zeros(Float64, N)
        Q_inward = zeros(Float64, N)
        P_inward[ctp:end], Q_inward[ctp:end], imin = schroed_inward_adams(l, E, V[ctp:end], mesh.r[ctp:end], mesh.rp[ctp:end])
        if imin > 1
            # The inward integration to reach the turning point, diverged
            @warn "Inward integration to reach the turning point diverged."
            break
        end

        # Normalize the inward to match the outward solution
        # P' = Q
        factor = P[ctp] / P_inward[ctp]
        if abs(factor) > 1e+9
            @warn "Normalization factor is too large"
            break
        end
        P_inward *= factor
        Q_inward *= factor

        # match the outward and inward solutions
        P[ctp+1:end] = P_inward[ctp+1:end]
        Q[ctp+1:end] = Q_inward[ctp+1:end]

        # E₂ ≈ E₁ + (Q₁⁻ - Q₁⁺) * P₁(ctp) / 2 ∫ P₁^2 dr
        S = integrate(P .^ 2, mesh.rp, method=:trapz7)
        dE = (Q[ctp] - Q_inward[ctp]) * P[ctp] / (2 * S)

        if abs(dE) < tol
            is_converged = true
            break
        end

        # dE = E₂ - E₁
        # dE < 0 means the new energy is lower, so the previous energy is too high
        # and is towards to the smaller energy
        if dE < 0
            Emax = E # towards to the smaller energy
        else
            Emin = E # towards to the larger energy
        end

        # if the dE prediction is out of the trust region, we still use bisection
        if E + dE < Emin || E + dE > Emax
            E = (Emax + Emin) / 2
            last_bisect = true
        else
            E += dE
            last_bisect = false
        end
    end

    # Normarlize the wave function
    S = integrate(P .^ 2, mesh.rp, method=:trapz7)
    S = sqrt(abs(S))

    if S > 0
        P /= S
        Q /= S
    else
        throw(ArgumentError("Wave function is too small, S = $S"))
    end

    E, P, Q, is_converged
end

function find_ctp(V::Vector{Float64}, E::Float64)::Int64
    N = length(V)
    for i in N:-1:1
        if V[i] < E
            return i
        end
    end
    N
end

function is_monotonic(P::Vector{Float64})::Bool
    N = length(P)
    for i in 2:N
        if abs(P[i]) < abs(P[i-1])
            return false
        end
    end
    return true
end

function count_nodes(P::Vector{Float64})::Int64
    N = length(P)
    nnodes = 0
    for i in 2:N
        # exclude the case P[i] == 0.0 once, to avoid double counting
        if sign(P[i]) != sign(P[i-1]) && P[i] != 0.0
            nnodes += 1
        end
    end

    nnodes
end

function get_min_idx(y::Vector{Float64})::Int64
    k = length(y)
    while abs(y[k-1]) < abs(y[k])
        k -= 1
        if k == 1
            break
        end
    end
    k = k - 1
    k
end