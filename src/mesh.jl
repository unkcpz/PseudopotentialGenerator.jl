using Dierckx: Spline1D
using Polynomials: Polynomial, fit, coeffs
import Polynomials: derivative

"""
    Mesh

A mesh object that stores the radial grid points and their derivatives.
To create a mesh object, use the constructor `Mesh(r_min, r_max, a, N)`.

The radial grid points are generated using the formula:

    r[i] = r_min + alpha * (exp(beta * (i - 1)) - 1)

    alpha = (r_max - r_min) / (exp(beta * N) - 1) 
    
    beta = log(a) / (N - 1).

The actual number of points in the mesh is N+1.
"""
struct Mesh
    r_min::Float64
    r_max::Float64
    a::Float64
    _N::Int64
    size::Int64
    r::Vector{Float64}
    rp::Vector{Float64}
    function Mesh(r_min::Float64, r_max::Float64, a::Float64, N::Int64)
        r = zeros(Float64, N+1)
        rp = zeros(Float64, N+1)

        if a < 0
            throw(ArgumentError("require a > 0"))
        elseif (abs(a - 1) < eps(typeof(a)))
            # Uniform grid if a = 1
            r = range(r_min, r_max, length=N+1)
            rp = fill((r_max - r_min) / N, N+1)
        else
            # Exponential grid
            if (N > 1)
                beta = log(a) / (N - 1)
                alpha = (r_max - r_min) / (exp(beta * N) - 1)
                for i = 1:N+1
                    r[i] = r_min + alpha * (exp(beta * (i - 1)) - 1)
                    rp[i] = alpha * beta * exp(beta * (i - 1))
                end
            elseif N == 1
                r = [r_min, r_max]
            else
                throw(ArgumentError("require N > 0"))
            end
        end

        new(r_min, r_max, a, N, N+1, r, rp)
    end
end

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


function mesh_exp_deriv2(mesh::Mesh)::Vector{Float64}
    rpp = mesh_exp_deriv2(mesh.r_min, mesh.r_max, mesh.a, mesh._N)
    rpp
end


function mesh_exp_deriv2(r_min::Float64, r_max::Float64, a::Float64, N::Int64)::Vector{Float64}
    rpp = zeros(Float64, N+1)
    if a < 0
        throw(ArgumentError("require a > 0"))
    elseif (abs(a - 1) < eps(typeof(a)))
        # Uniform grid if a = 1, rpp = 0 for all i
        nothing
    else
        # Exponential grid
        if (N > 1)
            beta = log(a) / (N - 1)
            alpha = (r_max - r_min) / (exp(beta * N) - 1)
            for i = 1:N+1
                rpp[i] = alpha * beta^2 * exp(beta * (i - 1))
            end
        elseif N == 1
            rpp = [0, 0]
        else
            throw(ArgumentError("require N > 0"))
        end
    end

    rpp
end

# TODO: will spline interpolation be better?
"""
    derivative(f, mesh, idx)

Compute the derivative of the function `f` at the point `idx` using the mesh `mesh`.
By fit the polynomial to the function values at the points `idx-5:idx+5` and compute the derivative at the point `idx`.

Can approximate the derivative to ~ 1e-8.
"""
function dfdr(f::Vector{Float64}, mesh::Mesh, idx::Int64)::Float64
    if idx < 6 || idx > mesh.size - 6
        throw(ArgumentError("The index must be greater than 6 and less than N-6"))
    end

    # Fit the polynomial
    rs = mesh.r[idx-5:idx+5]
    fs = f[idx-5:idx+5]
    p = fit(rs, fs, 6) |> p -> round.(coeffs(p), digits=10) |> Polynomial

    # Compute the derivative
    dp = derivative(p)

    dp(mesh.r[idx])
end

"""
    d2fdr2(f, mesh, idx)

Compute the second derivative of the function `f` at the point `idx` using the mesh `mesh`.
"""
function d2fdr2(f::Vector{Float64}, mesh::Mesh, idx::Int64)::Float64
    if idx < 11 || idx > mesh.size - 11
        throw(ArgumentError("The index must be greater than 6 and less than N-6"))
    end

    # Fit the polynomial
    rs = mesh.r[idx-5:idx+5]
    fs = f[idx-5:idx+5]
    p = fit(rs, fs, 6) |> p -> round.(coeffs(p), digits=12) |> Polynomial

    # Compute the derivative
    dp = derivative(p, 2)

    dp(mesh.r[idx])
end