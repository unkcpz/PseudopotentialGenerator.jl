function mesh_exp(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    mesh = zeros(Float64, N+1)
    @ccall libDFTATOM.mesh_exp(r_min::Ref{Float64}, r_max::Ref{Float64}, a::Ref{Float64}, N::Ref{Cint}, mesh::Ptr{Float64})::Vector{Float64}
    mesh
end

function mesh_exp_deriv(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    rp = zeros(Float64, N+1)
    @ccall libDFTATOM.mesh_exp_deriv(r_min::Ref{Float64}, r_max::Ref{Float64}, a::Ref{Float64}, N::Ref{Cint}, rp::Ptr{Float64})::Vector{Float64}
    rp
end

function mesh_exp_deriv2(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    rpp = zeros(Float64, N+1)
    @ccall libDFTATOM.mesh_exp_deriv2(r_min::Ref{Float64}, r_max::Ref{Float64}, a::Ref{Float64}, N::Ref{Cint}, rpp::Ptr{Float64})::Vector{Float64}
    rpp
end