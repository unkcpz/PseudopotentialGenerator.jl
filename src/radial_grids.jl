using Dierckx

export RGrid, primitive, gradient, hessian

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

function find_idx(R::RGrid, r::Float64)
    if r < R.r[1] || r > R.r[end]
        error("r is out of range")
    end
    idx = searchsortedlast(R.r, r)
    idx
end

function integrate(R::RGrid, f::Vector{Float64}, a::Float64, b::Float64)
    spl = Spline1D(R.r, f, k=5, bc="nearest")
    i = Dierckx.integrate(spl, a, b)
    i
end

function primitive(R::RGrid, f::Vector{Float64})
    i = integrate(R, f, R.r[1], R.r[end])
    i
end

function gradient(R::RGrid, f::Vector{Float64})
    # Using interpolation to calculate the derivative
    fd = Spline1D(R.r, f, k=5, bc="nearest")
    fd1 = Dierckx.derivative(fd, R.r)
    fd1
end

function hessian(R::RGrid, f::Vector{Float64})
    # Using interpolation to calculate the derivative
    grad_f = gradient(R, f)
    fd = Spline1D(R.r, grad_f, k=5, bc="nearest")
    fd2 = Dierckx.derivative(fd, R.r)
    fd2
end