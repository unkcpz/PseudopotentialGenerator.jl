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