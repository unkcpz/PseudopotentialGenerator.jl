"""
    lschfb solve the radial eigenvalue problem for a given potential.
        wrapping and calling lschfb from oncvpsp.
"""
LIBONCV_PATH = joinpath(@__DIR__, "../deps/liboncv.so")

function lschfb(n::Int, l::Int, ε::Float64, rgrid::Vector{Float64}, vloc::Vector{Float64}, Z::Int, srel::Bool)
    mmax = length(rgrid)
    nn = n
    ll = l
    ierr = Ref{Int32}(-99)  # for undefined
    ee = Ref{Float64}(ε)
    rr = rgrid
    vv = vloc
    uu = zeros(Float64, mmax)
    up = zeros(Float64, mmax)
    zz = Z
    mch = Ref{Int32}(0)
    srel = convert(Int, srel)
    
    @assert mmax == length(vloc)
    # lschfb use integer, logical kind=4 and real(kind(1.0d0))
    ccall(
        (:lschfb_, LIBONCV_PATH), 
        Nothing, 
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,srel, # same label and order as original function.  
    )
    uu, up, ee[], ierr[], mch[]
end

"""
    solve raidal Schrodinger equation of n, l orbital 
"""
function sol_orb(n::Int, l::Int, rgrid::Vector{Float64}, vloc::Vector{Float64}, Z::Int, srel::Bool; e_guess::Float64=0.0)
    uu, up, ee, ierr, mch = lschfb(n, l, e_guess, rgrid, vloc, Z, srel)
    if ierr != 0
        throw(ErrorException("lschfb not finished okay."))
    end

    uu, up, ee
end

# Thomas Fermi potential implemented from ONCVPSP src/tfapot.f90
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
function tf(Z::Int, r::Float64)
    zz = convert(Float64, Z)
    # TODO: this one is simple to implemented
    pot = ccall((:tfapot_, LIBONCV_PATH), Float64, (Ref{Float64}, Ref{Float64}), r, zz)
    pot
end

function tf(Z::Int, rgrid::Vector{Float64})
    pot = zeros(Float64, length(rgrid))
    for i in 1:length(rgrid)
        pot[i] = tf(Z, rgrid[i])
    end

    pot
end

# Copy the hartree.f90 is from PStudio code by D. Ceresoli which modified from vout.f90 of ONCVPSP
function hartree(Z::Int, rgrid::Vector{Float64}, ρ::Vector{Float64})
    rho = ρ
    gsize = length(ρ)
    mmax = Ref{Int32}(gsize)
    vo = zeros(Float64, gsize)
    zion = Float64(Z)
    rr = rgrid
    ccall((:hartree_, LIBONCV_PATH), Nothing, 
    (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}),
    rho, vo, zion, rr, mmax)

    vo
end