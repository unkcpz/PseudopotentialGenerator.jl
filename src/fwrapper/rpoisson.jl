function rpoisson_outward_pc(ρ::Vector{Float64}, r::Vector{Float64}, rp::Vector{Float64})::Vector{Float64}
    N = length(r)
    Vh = zeros(Float64, N)

    @ccall libDFTATOM.rpoisson_outward_pc(r::Ref{Float64}, rp::Ref{Float64}, ρ::Ref{Float64}, Vh::Ptr{Float64}, N::Ref{Int32})::Cvoid
    Vh
end