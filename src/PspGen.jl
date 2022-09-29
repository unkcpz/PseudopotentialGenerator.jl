module PspGen

include("./oncvpsp.jl")
include("./potential.jl")

function compute_rgrid(Z::Int; rmin::Float64=0.0001, rmax::Float64=10.0, N::Int=2000)
    xmin = log(Z * rmin)
    xmax = log(Z * rmax)
    xgrid = range(xmin, stop=xmax, length=N)
    dx = convert(Float64, xgrid.step)
    rgrid = @. exp(xgrid) / Z
    dr = @. dx * rgrid

    rgrid, dr
end

end