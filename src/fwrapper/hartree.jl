
function hartree(q::Float64, rgrid::Vector{Float64}, ρ::Vector{Float64})
    rho = ρ
    gsize = length(ρ)
    mmax = Ref{Int32}(gsize)
    vo = zeros(Float64, gsize)
    charge = Float64(q)
    rr = rgrid
    ccall((:hartree_, LIBONCV_PATH), Nothing, 
    (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}),
    rho, vo, charge, rr, mmax)

    vo
end