include("./oncvpsp.jl")

# Thomas Fermi potential implemented from ONCVPSP src/tfaport.f90
# ! generalized Thomas-Fermi atomic potential
# (Copy from the comment from ONCVPSP)
# !...to an article of N. H. March ( "The Thomas-Fermi Approximation in 
# ! Quantum Mechanics", Adv. Phys. 6, 1 (1957)). He has the formula,
# ! but it is not the result of his work. The original publication is: 
# !     R. Latter, Phys. Rev. 99, 510 (1955).
# ! He says that it''s an analytic fit to an improved calculation of the 
# ! potential distribution of a Thomas-Fermi atom without exchange first 
# ! performed by Miranda (C. Miranda, Mem. Acc. Italia 5, 285 (1934)).
# !                                 Alexander Seidl, TU Munich
function tf(Zion::Int, r::Float64)
    bb = (0.69395656/Zion)^(1/3)
    xx = r / bb
    xs = √xx

    tt = z/(1.0+xs*(0.02747 - xx*(0.1486 - 0.007298*xx)) + xx*(1.243 + xx*(0.2302 + 0.006944*xx)))

    if tt < 1.0
        tt = 1.0
    end

    pot = -tt / rr
    pot
end

function tf(Zion::Int, rgrid::Vector{Float64})
    pot = zeros(Float64, length(rgrid))
    for i in 1:length(rgrid)
        pot[i] = tf(Zion, rgrid[i])
    end

    pot
end

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
        # solve sch for every orbital
        
        # compute rho and mixing a new rho
        ρ_new = compute_rho(orbs, gsize)
        ρ_new = β * ρ_new + (1.0 - β) * ρ

        # compute hartree and xc potential


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
