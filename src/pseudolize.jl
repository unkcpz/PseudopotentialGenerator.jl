using Polynomials

function pseudolize(
    ϕ::Vector{Float64}, 
    nl::NamedTuple{(:n, :l), Tuple{Int64, Int64}}, 
    rc::Float64, 
    εl::Float64,
    vae::Vector{Float64},
    mesh::Mesh;
    method=:TM)
    if method == :TM
        v_pspot, ϕ_ps = pseudolize_TM(ϕ, nl, rc, εl, vae, mesh)
    else
        throw(ArgumentError("Unsupported pseudolization method $method"))
    end

    v_pspot, ϕ_ps
end

function pseudolize_TM(
    ϕ::Vector{Float64}, 
    nl::NamedTuple{(:n, :l), Tuple{Int64, Int64}}, 
    rc::Float64, 
    εl::Float64,
    vae::Vector{Float64},
    mesh::Mesh)

    # TODO: should add a checker for rc, if it is too large or too small
    # the NL solver may not converge

    # rc is user defined cutoff radius, it may not be in the grid
    # round it to the nearest grid point ≤ rc
    if rc > mesh.r[end] || rc < mesh.r[1]
        throw(ArgumentError("rc is out of the grid range"))
    end
    ic = findfirst(r -> r > rc, mesh.r) - 1
    rc = mesh.r[ic]
    # TODO: log it and use more detailt message
    println("New rc is: ", rc)

    # copy the ae wavefunction to ps wavefunction
    rϕ = mesh.r .* ϕ

    # make sure the ps wavefunction at rc is positive
    sign = 1.0
    if rϕ[ic] < 0.0
        rϕ .= -rϕ
        sign = -1.0
    end

    rϕ_ps = copy(rϕ)

    # Root finding for residual TM using NL solver
    # c2, c4, c6, c8, c10: coefficients that solved in linear equation
    cs = zeros(Float64, 5)

    function fixpoint_tm(c2, params)
        # static params
        ic = params.ic
        rc = params.rc
        vrc = params.vrc
        vprc = params.vprc
        vpprc = params.vpprc
        ae_norm = params.ae_norm
        εl = params.εl

        # The last NL condition (29g)
        c4 = - c2 ^ 2 / (2 * nl.l + 5)

        # solve the linear equation to get cs = [c0, c6, c8, c10, c12]
        # equations (29b), (29c), (29d), (29e), (29f) of TM paper PhysRevB.43.1993
        # rhs(1) = log(prc/rc**(l+1))
        # rhs(2) =   prc'/prc - (l + 1)/rc
        # rhs(3) =   2*(vrc - εl) - 2*(l + 1)/rc*rhs(2) - rhs(2)**2
        # rhs(4) =   2*vprc + 2*(l + 1)/rc**2*rhs(2) - 2*(l + 1)/rc*rhs(3) - 2*rhs(2)*rhs(3)
        # rhs(5) =   2*vpprc - 4*(l + 1)/rc**3*rhs(2) + 4*(l + 1)/rc**2*rhs(3) - 2*(l + 1)/rc*rhs(4) - 2*rhs(3)**2 - 2*rhs(2)*rhs(4)
        prc = rϕ[ic]
        dprc = dfdr(rϕ, mesh, ic)     # TODO: or f + f' * h?
        rhs = zeros(Float64, 5)
        rhs[1] = log(prc / rc ^ (nl.l + 1))
        rhs[2] = dprc / prc - (nl.l + 1) / rc
        rhs[3] = 2 * (vrc - εl) - 
                    2 * (nl.l + 1) / rc * rhs[2] - rhs[2] ^ 2
        rhs[4] = 2 * vprc + 2 * (nl.l + 1) / rc ^ 2 * rhs[2] - 
                    2 * (nl.l + 1) / rc * rhs[3] - 2 * rhs[2] * rhs[3]
        rhs[5] = 2 * vpprc - 4 * (nl.l + 1) / rc ^ 3 * rhs[2] + 
                    4 * (nl.l + 1) / rc ^ 2 * rhs[3] - 2 * (nl.l + 1) / rc * rhs[4] - 
                    2 * rhs[3] ^ 2 - 2 * rhs[2] * rhs[4]

        M_lhs = [ 1 rc^6  rc^8    rc^10     rc^12;
                  0 6rc^5  8rc^7   10rc^9    12rc^11;
                  0 30rc^4 56rc^6  90rc^8    132rc^10;
                  0 120rc^3 336rc^5 720rc^7   1320rc^9;
                  0 720rc^2 3360rc^4 10080rc^6 23760rc^8]

        # c2, c4
        rhs[1] -= c2 * rc^2   + c4 * rc^4
        rhs[2] -= c2 * 2 * rc + c4 * 4 * rc^3
        rhs[3] -= c2 * 2      + c4 * 12 * rc^2
        rhs[4] -= c2 * 0      + c4 * 24 * rc
        rhs[5] -= c2 * 0      + c4 * 24

        cs .= M_lhs \ rhs
        c0 = cs[1]
        c6 = cs[2]
        c8 = cs[3]
        c10 = cs[4]
        c12 = cs[5]

        # calculate rϕ_ps after linear solver
        p(x) = c0 + c2 * x^2 + c4 * x^4 + c6 * x^6 + c8 * x^8 + c10 * x^10 + c12 * x^12
        @. rϕ_ps[1:ic] = mesh.r[1:ic] ^ (nl.l+1) * exp(p(mesh.r[1:ic]))

        fx = integrate(rϕ_ps[1:ic] .^ 2, mesh.r[1:ic]) - ae_norm  #   (29a)

        # TODO: log it and use more detailt message
        println("fx: ", fx)

        fx
    end

    vrc = vae[ic]
    vprc = dfdr(vae, mesh, ic)
    vpprc = d2fdr2(vae, mesh, ic)
    ae_norm = integrate(rϕ[1:ic] .^ 2, mesh.r[1:ic])
    println("rϕ (ic): ", rϕ[ic])

    params = (; ic=ic, rc=rc, vrc=vrc, vprc=vprc, vpprc=vpprc, ae_norm=ae_norm, εl=εl)

    # The initial guess for c0, c12: coefficients that solved by nonliear solver
    cx0 = 0.0

    prob = NonlinearProblem(fixpoint_tm, cx0, params)
    sol = solve(prob, alg=Broyden(); abstol=1e-10, maxiters=200)
    c2 = sol.u
    c4 = - c2 ^ 2 / (2 * nl.l + 5)
    c0, c6, c8, c10, c12 = cs
    p = Polynomial([c0, 0, c2, 0, c4, 0, c6, 0, c8, 0, c10, 0, c12], :x)
    dpdx = derivative(p)
    d2pdx2 = derivative(p, 2)

    # Get the ps wavefunction and recover the sign to AE wavefunction
    ϕ_ps = rϕ_ps ./ mesh.r
    ϕ_ps .= sign * ϕ_ps

    # using p(x) to represent to avoid d2fdr2 with interpolation
    v_pspot = copy(vae)
    @. v_pspot[1:ic] = εl + 
            (nl.l + 1) * dpdx(mesh.r[1:ic]) / mesh.r[1:ic] + 
            (d2pdx2(mesh.r[1:ic]) + dpdx(mesh.r[1:ic]) ^ 2) / 2

    v_pspot, ϕ_ps
end
