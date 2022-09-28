"""
    lschfb solve the radial eigenvalue problem for a given potential.
        wrapping and calling lschfb from oncvpsp.
"""
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
        (:lschfb_, joinpath(@__DIR__, "../deps/liboncv.so")), 
        Nothing, 
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        nn,ll,ierr,ee,rr,vv,uu,up,zz,mmax,mch,srel, # same label and order as original function.  
    )
    uu, up, ee[], ierr[], mch[]
end

"""
    solve raidal Schrodinger equation of n, l orbital 
"""
function sol_orb(n::Int, l::Int, rgrid::Vector{Float64}, vloc::Vector{Float64}, Z::Int, srel::Bool)
    e_guess = 0.0
    uu, up, ee, ierr, mch = lschfb(n, l, e_guess, rgrid, vloc, Z, srel)
    if ierr != 0
        throw(ErrorException("lschfb not finished okay."))
    end

    uu, up, ee
end

function tfapot(rgrid::Vector{Float64}, Z::Int)

end