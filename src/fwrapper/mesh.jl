"""
    mesh(r_min, r_max, a, N) -> mesh

Create a mesh with `npts` points between `r_min` and `r_max` using an exponential grid with parameter `a`.
"""
function mesh(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    mesh = zeros(Float64, N+1)
    @ccall libDFTATOM.__wrapper_MOD_mesh_exp(r_min::Ref{Float64}, r_max::Ref{Float64}, a::Ref{Float64}, N::Ref{Int32}, mesh::Ptr{Float64})::Vector{Float64}
    mesh
end