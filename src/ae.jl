include("./oncvpsp.jl")

# Running atomic DFT calcualtion
function scf(Zion::Int, rgrid::Vector{Float64}; β::Float64=0.25, tol::Float64=1e-8)
    # v_ks = v_h + v_xc + v_ion
    # initial v_ks to thomas fermi atom potential
    v_ks = tf(Zion, rgrid)
    gsize = length(rgrid)

    ρ = zeros(Float64, gsize)
    E = 0.0

    i = 1
    flag_conv = false
    while i < max_iter
        ## solve sch for every orbital
        
        # compute rho and mixing a new rho
        ρ_new = compute_rho(orbs, gsize)
        ρ_new = β * ρ_new + (1.0 - β) * ρ

        # compute hartree and xc potential
        ## hartree


        # compute energy terms
        E_new = E_k + E_h + E_xc + E_vxc + E_ion
        δ = E_new - E
        E = E_new

        # Print the log

        # check convergence on energy
        if δ < tol
            flag_conv = true
            break
        end
    end

    flag_conv, ρ, E
end

"""
ρ = ∑ f * <ψ|ψ> in 1D radial 
"""
function compute_rho(orbs::Vector{Orbital}, gsize::Int)
    ρ = zeros(Float64, gsize)
    for orb in orbs
        ρ += orb.occ * orb.ur * orb.ur
    end

    ρ
end

# give hartree potential from density
function compute_hartree(ρ::Vector{Float64})
end