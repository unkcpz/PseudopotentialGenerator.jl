function compute_ld(l::Int, Z::Int, ε::Float64, mesh::Mesh, vae::Vector{Float64}, rc::Float64)
    # rc is user defined cutoff radius, it may not be in the grid
    # round it to the nearest grid point ≤ rc
    if rc > mesh.r[end] || rc < mesh.r[1]
        throw(ArgumentError("rc is out of the grid range"))
    end
    ic = findfirst(r -> r > rc, mesh.r) - 1
    rc = mesh.r[ic]

    # TODO: make the final integrat point setable to save time
    P, Q, _ = schroed_outward_adams(l, Z, ε, vae, mesh.r, mesh.rp)

    # Normarlize the wave function
    S = integrate(P .^ 2, mesh.rp, method=:trapz7)
    S = sqrt(abs(S))

    if S > 0
        P /= S
        Q /= S
    else
        throw(ArgumentError("Wave function is too small, S = $S"))
    end

    # P = rR(r); Q = rR'(r) + R(r)
    # logder = R'(r) / R(r) = Q / P - 1 / r
    return Q[ic] / P[ic] - 1 / mesh.r[ic]
end

function compute_atanld(l::Int, Z::Int, mesh::Mesh, vae::Vector{Float64}, rc::Float64; window=[-20.0, 20.0], δ=0.01)
    x = window[1]:δ:window[2]
    y = [atan(compute_ld(l, Z, ε, mesh, vae, rc)) for ε in x]

    # y can have a π shift since arctan is not unique, we need to find the correct shift
    for idx in eachindex(y)[2:end]
        s = y[idx] - y[idx - 1]
        s_factor = round(s / π)
        y[idx] -= s_factor * π
    end

    # For clarity, I shift y to y + l * π
    y .+= l * π

    return x, y
end