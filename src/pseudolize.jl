using Polynomials: derivative

"""
    pseudolize(ϕ, ε, vae, mesh, rc; method=:TM)

Pseudolize on the pseudo configurations. Loop over the orbitals and do the pseudolization.
The function also manage to unscreened the pseudo potential and convert it to the 
Kleinman-Bylander form.

Methods supported:
    - `:TM`: Troullier-Martins pseudopotential
    - `:TM22`: (TODO) Troullier-Martins pseudopotential with 22 coefficients (doi:10.1088/1361-648X/aac85d)
    - `:RRKJ`: (TODO) Rappe-Rabe-Kaxiras-Joannopoulos pseudopotential (https://doi.org/10.1103/PhysRevB.41.1227)
    - `:BHS`: (TODO) Bachelet-Hamann-Schlüter pseudopotential (doi:10.1103/PhysRevB.26.4199)
    - more...
"""
function pseudolize(
    ae_info,    # TODO: type it
    mesh::Mesh, 
    rc::Dict{NamedTuple{(:n, :l), Tuple{Int64, Int64}}, Float64}; 
    method=:TM,
    kbform=true)::Potential

    # To store the results from pseudolization for post processing
    v_pspot_dict = Dict{NamedTuple{(:n, :l), Tuple{Int64, Int64}}, Vector{Float64}}()
    ϕ_ps_dict = Dict{NamedTuple{(:n, :l), Tuple{Int64, Int64}}, Vector{Float64}}()

    # loop over rc which defined the target orbitals to pseudolize
    # The ϕs and εs are the all wavefunctions and eigenvalues from AE calculation
    for (nl, _rc) in rc
        ϕ = ae_info.ϕs[nl]
        ε = ae_info.ε_lst[nl]
        vae = ae_info.vae
        v_pspot, ϕ_ps = pseudolize(nl, ϕ, ε, vae, mesh, _rc; method=method)
        v_pspot_dict[nl] = v_pspot

        # normalize the pseudo wavefunction and store it
        ϕ_ps = ϕ_ps / sqrt(integrate(ϕ_ps .^ 2, mesh.r))
        ϕ_ps_dict[nl] = ϕ_ps
    end

    # Unscreened the pseudo by subtracting vxc[ρ] and vh[ρ] from each v_pspot.
    # ρ is ρ_ps charge density from the pseudo wavefunctions
    # The pseudo-wave of every orbital is get from AE wave by "fitting" so for each of them the vxc and vh are counted.
    # Note for myself to understanding, I was thinking substructing vxc and vh from every l is dupulicated, but it is correct. Consider if the v_pspot is the same as v_ae, then the potential will be counting multiple times. 

    # compute the density from pseudo wavefunctions
    ρ = zeros(Float64, length(mesh.r))
    for (nl, ϕ_ps) in ϕ_ps_dict
        @. ρ += ϕ_ps ^ 2 * ae_info.occs[nl]
    end

    # compute the vxc and vh from the pseudo density
    vxc = compute_vxc(ρ, ae_info.xc)
    vh = compute_vh(ρ, mesh)

    # unscreen the pseudo potential
    for k in keys(v_pspot_dict)
        @. v_pspot_dict[k] -= vxc.v_xc + vh
    end    

    v_sl = SemiLocalPotential(v_pspot_dict)

    if kbform
        v_kb = sl2kb(v_sl)
        return v_kb
    else
        return v_sl
    end
end

function sl2kb(v_sl::SemiLocalPotential)::KBFormPotential
    # Convert the semi-local potential to Kleinman-Bylander form

    # The local potential can, in principle, be arbitrarily chosen, but since the summation of nonlocal part need to be truncated at some value of l, 
    # the local potential should be chosen such that it adequately reproduces
    #the atomic scattering for all the higher angularmomentum channels.
    # TODO: more versatile way to choose the local potential, check https://arxiv.org/pdf/1811.00959.pdf see anything can be revived from the paper.
    # use l=0 for the local potential for now
    v_pspot_local = zeros(Float64, length(v_sl.v))
    for (nl, v_pspot) in v_pspot_dict
        if nl.l == 0
            v_pspot_local .= v_pspot
            break
        end
    end

    # XXX ekb, proj_kb are not implemented yet
end

"""
    pseudolize(nl, ϕ, ε, vae, mesh, rc; method=:TM)

Pseudolize on the given orbital.
"""
function pseudolize(
    nl::NamedTuple{(:n, :l), Tuple{Int64, Int64}}, 
    ϕ::Vector{Float64}, 
    ε::Float64,
    vae::Vector{Float64},
    mesh::Mesh,
    rc::Float64; 
    method=:TM)
    if method == :TM
        v_pspot, ϕ_ps = pseudolize_TM(nl, ϕ, ε, vae, mesh, rc)
    else
        throw(ArgumentError("Unsupported pseudolization method $method"))
    end

    v_pspot, ϕ_ps
end

function pseudolize_TM(
    nl::NamedTuple{(:n, :l), Tuple{Int64, Int64}}, 
    ϕ::Vector{Float64}, 
    ε::Float64,
    vae::Vector{Float64},
    mesh::Mesh,
    rc::Float64)

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

    # ϕ is the radial wavefunction R(r) in most paper representation
    rϕ = mesh.r .* ϕ

    # make sure the ps wavefunction at rc is positive
    sign = 1.0
    if rϕ[ic] < 0.0
        rϕ .= -rϕ
        sign = -1.0
    end

    # copy the ae wavefunction to ps wavefunction
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
        ε = params.ε

        # The last NL condition (29g)
        c4 = - c2 ^ 2 / (2 * nl.l + 5)

        # solve the linear equation to get cs = [c0, c6, c8, c10, c12]
        # equations (29b), (29c), (29d), (29e), (29f) of TM paper PhysRevB.43.1993
        # rhs(1) = log(prc/rc**(l+1))
        # rhs(2) =   prc'/prc - (l + 1)/rc
        # rhs(3) =   2*(vrc - ε) - 2*(l + 1)/rc*rhs(2) - rhs(2)**2
        # rhs(4) =   2*vprc + 2*(l + 1)/rc**2*rhs(2) - 2*(l + 1)/rc*rhs(3) - 2*rhs(2)*rhs(3)
        # rhs(5) =   2*vpprc - 4*(l + 1)/rc**3*rhs(2) + 4*(l + 1)/rc**2*rhs(3) - 2*(l + 1)/rc*rhs(4) - 2*rhs(3)**2 - 2*rhs(2)*rhs(4)
        prc = rϕ[ic]
        dprc = dfdr(rϕ, mesh, ic)     # TODO: or f + f' * h?
        rhs = zeros(Float64, 5)
        rhs[1] = log(prc / rc ^ (nl.l + 1))
        rhs[2] = dprc / prc - (nl.l + 1) / rc
        rhs[3] = 2 * (vrc - ε) - 
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

    params = (; ic=ic, rc=rc, vrc=vrc, vprc=vprc, vpprc=vpprc, ae_norm=ae_norm, ε=ε)

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
    @. v_pspot[1:ic] = ε + 
            (nl.l + 1) * dpdx(mesh.r[1:ic]) / mesh.r[1:ic] + 
            (d2pdx2(mesh.r[1:ic]) + dpdx(mesh.r[1:ic]) ^ 2) / 2

    v_pspot, ϕ_ps
end
