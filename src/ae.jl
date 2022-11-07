using Libxc

include("./oncvpsp.jl")

mutable struct Orbital
    n::Int
    l::Int
    occ::Float64
    Z::Int
    ε::Float64
    rgrid::Vector{Float64}
    ur::Vector{Float64}
    function Orbital(n, l, occ, Z, rgrid)
        ε = -0.5 * Z^2 / n^2
        ur = zeros(Float64, length(rgrid))
        new(n, l, occ, Z, ε, rgrid, ur)
    end
end

# TODO configuration parsed to orbitals

# Running atomic DFT calcualtion
function scf!(
    Z::Int, 
    orbs::Vector{Orbital}, 
    rgrid::Vector{Float64}, 
    dr::Vector{Float64},
    xc::Tuple{Symbol, Symbol},
    srel::Bool;
    max_iter::Int=300, β::Float64=0.25, tol::Float64=1e-8
)
    # v_ks = v_h + v_xc + v_ion
    # initial v_ks to thomas fermi atom potential
    v_ks = tf(Z, rgrid)
    gsize = length(rgrid)

    rho = zeros(Float64, gsize)
    rho_old = zeros(Float64, gsize)
    E_tot = 0.0
    E_old = 0.0

    i = 1
    flag_conv = false
    while i < max_iter
        ## solve sch for every orbital
        sol!(orbs, rgrid, dr, v_ks, Z, srel)
        
        # compute rho and mixing a new rho
        rho = compute_rho(orbs, rgrid)
        if i > 1
            rho = β * rho + (1.0 - β) * rho_old
        end
        
        rho_old .= rho # copy and store for next iteration

        # comptue kenitic energy
        E_k = compute_Ek(orbs)

        # compute hartree and xc potential
        ## hartree
        E_h, v_h = compute_hartree(rho, Z, rgrid, dr)

        ## xc
        # functional = Functional(:lda_x, n_spin=1) # XXX should be a list and evaluate x and c respectively
        E_xc, v_xc, E_vxc = compute_xc(rho, rgrid, dr, xc)

        # ion
        E_ion, v_ion = compute_ion(Z, rho, rgrid, dr)

        # update v_ks
        v_ks = v_ion + v_h + v_xc

        # compute energy terms
        E_tot = E_k - E_h + E_xc - E_vxc
        E_k = E_tot - E_ion - E_h - E_xc
        δ = E_tot - E_old
        E_old = E_tot

        println("E_k: ", E_k)
        println("E_h: ", E_h) 
        println("E_xc: ", E_xc)
        println("E_vxc: ", E_vxc)
        println("E_ion: ", E_ion)

        # Print the log

        # update v_ks
        # v_ks = v

        # check convergence on energy
        if abs(δ) < tol
            flag_conv = true
            break
        end
    end

    flag_conv, rho, E_tot
end

# solve and update the orbitals
function sol!(orbs::Vector{Orbital}, rgrid::Vector{Float64}, dr::Vector{Float64}, v::Vector{Float64}, Z::Int, srel::Bool)
    for orb in orbs
        n, l, occ, e_guess = orb.n, orb.l, orb.occ, orb.ε
        if occ < 0.0    # unbound state
            orb.ur = zeros(Float64, size(rgrid))
            orb.ε = 0.0
            continue
        end

        try
            # XXX: check if e_guess from previous step is needed or cause issue?
            orb.ur, _, orb.ε = sol_orb(n, l, rgrid, v, Z, srel, e_guess=e_guess)
        catch exc
            if isa(exc, ErrorException)
                orb.ur = zeros(Float64, size(rgrid))
                orb.ε = 0.0
            end
        finally
            # normalize the wavefunction
            urur = itrap(orb.ur .* orb.ur, dr)
            @. orb.ur /= sqrt(urur)
        end


    end
end


"""
ρ = ∑ f * <ψ|ψ> in 1D radial 
"""
function compute_rho(orbs::Vector{Orbital}, rgrid::Vector{Float64})
    rho = zeros(Float64, length(rgrid))
    for orb in orbs
        @. rho += orb.occ * orb.ur * orb.ur
    end
    @. rho = rho / (4π * rgrid * rgrid)
    rho
end

function compute_Ek(orbs::Vector{Orbital})
    # sum of state energy (with occupation)
    E_k = 0.0
    for orb in orbs
        E_k += orb.occ * orb.ε
    end

    E_k
end

# The energy E = 0.5 * 4π ∫ ρ(r)v(r) r^2 dr in radial space
# give hartree potential from density
# return energy and potential
function compute_hartree(rho::Vector{Float64}, Z::Int, rgrid::Vector{Float64}, dr::Vector{Float64})
    v_h = hartree(Z, rgrid, 4π * rho)   # why 4π??
    E_h = 0.5 * integrate(rho, v_h, rgrid, dr)
    E_h, v_h
end

# return ion energy and potential
function compute_ion(Z::Int, rho::Vector{Float64}, rgrid::Vector{Float64}, dr::Vector{Float64})
    v_ion = @. -Z / rgrid
    E_ion = integrate(rho, v_ion, rgrid, dr)

    E_ion, v_ion
end

function compute_rgrid(Z::Int; rmin::Float64=0.00001, rmax::Float64=100.0, N::Int=2000)
    xmin = log(Z * rmin)
    xmax = log(Z * rmax)
    xgrid = range(xmin, stop=xmax, length=N)
    dx = convert(Float64, xgrid.step)
    rgrid = @. exp(xgrid) / Z
    dr = @. dx * rgrid

    rgrid, dr
end

# xc energy and potential
function compute_xc(rho::Vector{Float64}, rgrid, dr, xc::Tuple{Symbol, Symbol})
    # compute exchange
    res = evaluate(Functional(xc[1], n_spin=1), rho=rho)

    exc = res.zk
    vxc = vec(res.vrho)

    res = evaluate(Functional(xc[2], n_spin=1), rho=rho)
    exc += res.zk
    vxc += vec(res.vrho)

    E_xc = integrate(rho, exc, rgrid, dr)   # ?? why still times charge density?
    E_vxc = integrate(rho, vxc, rgrid, dr)  # integral of rho * vxc ??

    E_xc, vxc, E_vxc
end

"""
    integral over space. For energy etc.
"""
function integrate(rho, v, rgrid, dr)
    integrand = @. v * rho * rgrid ^ 2
    E = itrap(integrand, dr) * 4π

    E
end
    
# trapezoidal method for uniform grid
# dx is the interval
function itrap(f::Vector{T}, dx::T)::T where {T <: Real}
    a = dx * sum(i -> f[i], 1:length(f) - 1)
    b = dx/2 * (f[begin] + f[end])

    a + b
end

# trapezoidal method for non-uniform grid
# dr <- r <- f
function itrap(f::Vector{T}, dr::Vector{T})::T where {T<: Real}
    @assert length(f) == length(dr)
    res = sum(i -> f[i] * dr[i], 1:length(f))
    res -= (f[begin] * dr[begin] + f[end] * dr[end]) / 2

    res
end

function tf(Z::Int, r::Float64)
    b = (0.69395656/Z) ^ (1.0/3.0)
    x = r / b
    xs = √x

    t = Z / (1.0 + xs * (0.02747 - x*(0.1486 - 0.007298*x)) + x * (1.243 + x*(0.2302 + 0.006944*x)))

    if t < 1.0
        t = 1.0
    end

    -t / r
end

function tf(Z::Int, rgrid::Vector{Float64})
    pot = zeros(Float64, length(rgrid))
    for i in 1:length(rgrid)
        pot[i] = tf(Z, rgrid[i])
    end

    pot
end