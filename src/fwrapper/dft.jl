function get_Vxc(ρ, r::Vector{Float64})
    N = length(r)
    vxc = zeros(Float64, N)
    exc = zeros(Float64, N)
    relat = false

    c = SPEED_OF_LIGHT

    @ccall libDFTATOM.get_vxc(r::Ptr{Float64}, ρ::Ptr{Float64}, relat::Ref{Cint}, c::Ref{Float64}, exc::Ptr{Float64}, vxc::Ptr{Float64}, N::Ref{Int32})::Cvoid

    (; v_xc=vxc, ε_xc=exc)
end