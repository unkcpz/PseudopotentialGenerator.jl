export RGrid

struct RGrid
    r::Vector{Float64}
    dr::Vector{Float64}
    N::Int64
    RGrid(r::Vector{Float64}, dr::Vector{Float64}) = new(r, dr, length(r))
end

function RGrid(Zion::Int64; dx=8e-3, xmin=-8, rmax=100)
    xmax = log(Zion * rmax)
    N = (xmax - xmin) / dx + 1
    N = 2 * (N / 2) + 1

    r = exp.(xmin .+ dx .* (0:N-1)) ./ Zion
    dr = r .* dx
    RGrid(r, dr)
end