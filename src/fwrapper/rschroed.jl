function schroed_outward_adams(l::Int64, Z::Int64, E::Float64, V::Function, r::Vector{Float64}, rp::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    P = zeros(Float64, N)
    Q = zeros(Float64, N)
    Vs = V.(r)
    imax = @ccall libDFTATOM.schroed_outward_adams(l::Ref{Cint}, Z::Ref{Cint}, E::Ref{Float64}, r::Ptr{Float64}, rp::Ptr{Float64}, Vs::Ptr{Float64}, P::Ptr{Float64}, Q::Ptr{Float64}, N::Ref{Cint})::Int64

    P, Q, imax
end

function schroed_inward_adams(l::Int64, E::Float64, V::Function, r::Vector{Float64}, rp::Vector{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Int64}
    N = length(r)
    P = zeros(Float64, N)
    Q = zeros(Float64, N)
    Vs = V.(r)
    imin = @ccall libDFTATOM.schroed_inward_adams(l::Ref{Cint}, E::Ref{Float64}, r::Ptr{Float64}, rp::Ptr{Float64}, Vs::Ptr{Float64}, P::Ptr{Float64}, Q::Ptr{Float64}, N::Ref{Cint})::Int64

    P, Q, imin
end