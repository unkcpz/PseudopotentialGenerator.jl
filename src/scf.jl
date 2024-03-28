using NonlinearSolve
using NLsolve
using Libxc: Functional, evaluate

mutable struct SavedStep
    ε_lst::Vector{Float64}
    Ys::Matrix{Float64}
    E_tot::Float64
    v_tot::Vector{Float64}
    v_h::Vector{Float64}
    v_xc::Vector{Float64}
    ε_xc::Vector{Float64}
    ρ::Vector{Float64}

    function SavedStep(Nr, Norb)
        ε_lst = zeros(Float64, Norb)
        Ys = zeros(Float64, Norb, Nr)
        E_tot = 0.0
        v_tot = zeros(Float64, Nr)
        v_h = zeros(Float64, Nr)
        v_xc = zeros(Float64, Nr)
        ε_xc = zeros(Float64, Nr)
        ρ = zeros(Float64, Nr)

        new(ε_lst, Ys, E_tot, v_tot, v_h, v_xc, ε_xc, ρ)
    end
end

# define the SCF problem
function self_consistent_field(Z::Int64, mesh::Mesh, orbs::Vector{Orbital}; mixing_beta::Float64=0.4, abstol::Float64=1e-6, maxiters::Int64=100)
    # use nonliear solver to solve the problem to find fixed point of: 
    # f_ks => residual(V) = KS(V) - V

    # f is the residual of the potential (non-Coulombic part)
    # V_total = V_coulomb + V_noncoulomb

    # The initial guess for the potential: TF potential
    v_tot = @. thomas_fermi_potential(Z, mesh.r) 
    v_coulomb = @. coulomb_potential(Z, mesh.r)
    r_size = length(mesh.r)

    v0 = v_tot - v_coulomb
    
    # Parameters to define the atom status, which are parameters should not 
    # changing during the SCF iteration
    params = (Z=Z, mesh=mesh, orbs=orbs, v_coulomb=v_coulomb)

    # For storing the step information
    saved = SavedStep(r_size, length(orbs))

    # initial v_tot
    saved.v_tot = v_tot

    # Use the hydrogen-like atom energy as the initial guess
    ε_lst = [hydrogen_like_energy(orb_idx, Z) for orb_idx in eachindex(orbs)]
    saved.ε_lst = ε_lst
    iter = 0

    function fixpoint_ks(vin, params)
        println("iter = ", iter)
        println("vin = ", vin[50])
        println("v_tot = ", saved.v_tot[50])
        println("v_coulomb = ", params.v_coulomb[50])
        println("ks_energies = ", saved.ε_lst)
        println("v_h = ", saved.v_h[50])
        println("v_xc = ", saved.v_xc[50])
        println("rho = ", saved.ρ[50])

        ρ = v2ρ!(vin, params, saved)
        vout = ρ2v!(ρ, params, saved)

        println("vout = ", vout[50])

        saved.ρ = ρ

        # calculate the total energy
        # TODO: log the total energy
        #E_tot = total_energy(ρ, params, saved)
        #saved.E_tot = E_tot

        iter += 1

        Δv = vin - vout
        Δv
    end

    prob = NonlinearProblem(fixpoint_ks, v0, params)
    solve(prob, alg=NLsolveJL(; method=:anderson, beta=mixing_beta); abstol=abstol, maxiters=maxiters)

    # total energy
    (; E_tot=saved.E_tot, v_tot=saved.v_tot, ρ=saved.ρ)
end

function v2ρ!(v_in::Vector{Float64}, p, saved)::Vector{Float64}
    # given the potential, calculate the density by solving the KS equations
    # the eigenvalues are updated in place
    r = p.mesh.r
    N = length(r)
    rho = zeros(Float64, N)

    # loop over the orbitals and solve the radial KS equations
    for idx in eachindex(p.orbs)
        # solve the KS equation
        orb = p.orbs[idx]
        n, l, f = orb.n, orb.l, orb.f

        Emax = orb.emax
        Emin = orb.emin
        Eini = saved.ε_lst[idx]

        v_tot = v_in + p.v_coulomb
        ε, P, _, is_converged = solve_radial_eigenproblem(n, l, p.Z, v_tot, p.mesh; tol=1e-9, max_iter=200, E_window=[Emin, Emax], E_ini = Eini, rel = false, perturb = true)
        println("ε = ", ε)

        if !is_converged
            throw(ArgumentError("Failed to solve the radial eigenproblem"))
        end

        # update the eigenvalues
        saved.ε_lst[idx] = ε

        # only rel = false is supported
        Y = P ./ r

        # update rho and orbitals
        #println("f: ", f)
        @. rho += Y ^ 2 * f
        saved.Ys[idx, :] = Y
    end

    # normalize the density
    rho = rho / (4 * π)
    #println("rho = ", rho[10:100:500])
    rho
end

function ρ2v!(ρ::Vector{Float64}, p, saved)::Vector{Float64}
    # Compute and store the step information: v_xc, v_h
    v_h = compute_vh(ρ, p.mesh)
    saved.v_h = v_h

    res = compute_vxc(ρ, :lda)
    v_xc = res.v_xc
    ε_xc = res.ε_xc
    saved.v_xc = v_xc
    saved.ε_xc = ε_xc

    E_tot = total_energy(ρ, p, saved)
    saved.E_tot = E_tot

    v_out = v_xc + v_h
    saved.v_tot = v_out + p.v_coulomb
    v_out
end

function compute_vxc(ρ::Vector{Float64}, xc::Symbol)
    # compute the exchange-correlation potential
    if xc == :lda
        func_x = Functional(:lda_x)
        func_c = Functional(:lda_c_vwn)
    else
        throw(ArgumentError("Unsupported exchange-correlation functional $xc"))
    end

    result_x = evaluate(func_x, rho=ρ)
    result_c = evaluate(func_c, rho=ρ)

    v_xc = result_x.vrho[:] + result_c.vrho[:]
    ε_xc = result_x.zk + result_c.zk
    (; v_xc=v_xc, ε_xc=ε_xc)
end

function compute_vh(ρ, mesh::Mesh)::Vector{Float64}
    # compute the Hartree potential
    vh = rpoisson_outward_pc(ρ, mesh)
    vh
end

function total_energy(ρ, p, saved)::Float64
    # compute the total energy
    r = p.mesh.r
    rp = p.mesh.rp
    occupations = [orb.f for orb in p.orbs]
    ks_energies = saved.ε_lst
    v_tot = saved.v_tot
    v_h = saved.v_h
    ε_xc = saved.ε_xc

    ∑ks_eneries = sum(occupations .* ks_energies)
    T_s = ∑ks_eneries - 4π * integrate(ρ .* v_tot .* r .^ 2, rp)

    E_ee = 2π * integrate(ρ .* v_h .* r .^ 2, rp)
    E_en = 4π * integrate(ρ .* p.v_coulomb .* r .^ 2, rp)
    E_xc = 4π * integrate(ρ .* ε_xc .* r .^ 2, rp)

    E_tot = T_s + E_ee + E_en + E_xc
    E_tot
end