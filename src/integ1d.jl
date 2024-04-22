"""
    integrate(y, x; method::Symbol)

Integrate the function y(x) over x using the specified method (default is F_Trapezoidal7).
The method must be one of the following (f_ prefix means the implementation is in Fortran):
    * `trapz1`
    * `trapz3`
    * `trapz5`
    * `trapz7`
    * `simpson`
    * `adams`
"""
# TODO: FastGaussQuadrature.jl check https://discourse.julialang.org/t/numerical-integration-only-with-evaluated-point/63552/5
function integrate(y, rp; method::Symbol=:trapz7)::Float64
    N = length(rp)
    if N < 8
        throw(ArgumentError("Length of integrated f must be at least 8"))
    end

    if method == :trapz1
        s = trapz1(y, rp)
    elseif method == :trapz3
        s = trapz3(y, rp)
    elseif method == :trapz5
        s = trapz5(y, rp)
    elseif method == :trapz7
        s = trapz7(y, rp)
    elseif method == :simpson
        s = simpson(y, rp)
    elseif method == :adams
        s = adams(y, rp)
    else
        throw(ArgumentError("Invalid method: $method"))
    end

    s
end

function trapz1(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (g[1] + g[N]) / 2
    s += sum(g[2:N-1])
    s
end

function trapz3(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        9 * (g[1] + g[N]) + 
        28 * (g[2] + g[N-1]) + 
        23 * (g[3] + g[N-2])
    ) / 24
    s += sum(g[4:N-3])
    s
end

function trapz5(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        475 * (g[1] + g[N]) +
        1902 * (g[2] + g[N-1]) +
        1104 * (g[3] + g[N-2]) +
        1586 * (g[4] + g[N-3]) +
        1413 * (g[5] + g[N-4])
    ) / 1440
    s += sum(g[6:N-5])
    s
end

function trapz7(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = (
        36799 * (g[1] + g[N]) +
        176648 * (g[2] + g[N-1]) +
        54851 * (g[3] + g[N-2]) +
        177984 * (g[4] + g[N-3]) +
        89437 * (g[5] + g[N-4]) +
        130936 * (g[6] + g[N-5]) +
        119585 * (g[7] + g[N-6])
    ) / 120960
    s += sum(g[8:N-7])
    s
end

function simpson(rp::Vector{Float64}, y::Vector{Float64})::Float64
    N = length(rp)
    g = y .* rp
    s = 0
    for i in 2:2:N-1
        s += g[i-1] + 4 * g[i] + g[i+1]
    end
    s /= 3
    if N % 2 == 0
        # If N is even, add the last slice
        s += (5 * g[N] + 8 * g[N-1] - g[N-2]) / 12
    end
    s
end
