# This file include the implementation copy from PGEN
# just for legacy test
using Dierckx: Spline1D

"""
    midpoints(fs, r)

Get the midpoints of the function values fs at the points r / 2, using interpolate.
"""
function midpoints(fs::Vector{Float64}, r::Vector{Float64})::Vector{Float64}
    if length(fs) != length(r)
        throw(ArgumentError("You must provide the same number of elements in fs and r"))
    end

    # Interpolate the function
    itp = Spline1D(r, fs)

    fs_mid = zeros(Float64, length(r)-1)

    for i = 1:length(r)-1
        mid = (r[i] + r[i+1]) / 2
        fs_mid[i] = itp(mid)
    end

    fs_mid
end

function midpoints(r::Vector{Float64})::Vector{Float64}
    rmid = zeros(Float64, length(r)-1)
    for i = 1:length(r)-1
        rmid[i] = (r[i] + r[i+1]) / 2
    end
    rmid
end

function integrate(y, rp; method::Symbol=:trapz7)::Float64
    N = length(rp)
    if N < 8
        throw(ArgumentError("Length of integrated f must be at least 8"))
    end

    if method == :trapz7
        s = trapz7(y, rp)
    else
        throw(ArgumentError("Invalid method: $method"))
    end

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