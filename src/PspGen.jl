module PspGen

# struct RadialGrid
#     Z::Int  # atom number zeta
#     N::Int  # number of points
#     rmin::Float64   # min of radial
#     rmax::Float64   # max of radial
#     xmin::Float64   # min of x in log mesh; xmin = log(Z * rmin)
#     xmax::Float64   # max of x in log mesh; xmax = log(Z * rmax)
#     x
#     dx
#     r
#     dr
# end

# struct Orbital
#     n::Int  # principal quantum number
#     l::Int  # angular quantum number
#     f::Float64  # occupation; f<0.0 for unbound states
#     rg::RadialGrid    # 
#     ϕr::Vector{Float64}  # wavefunction ϕ = ϕr / r
#     ev::Float64 # eigenvalue
#     function Orbital(n, l, rg, potential, reltype)
#         ϕr, ev = compute_orbital(n, l, rg.r, pot, Z, srel)
#         new(n, l, f, rg, ϕr, ev)
#     end
# end

# """
# Finds bound states of an all-electron atomic potential using
# Pauli-type scalar-relativistic Schroedinger equation
# https://arxiv.org/pdf/1207.5752.pdf
# """
# function compute_orbital(n, l, r, pot, Z, srel)

struct QuantumNumber
    n::Int  # principal quantum number
    l::Int  # angular quantum number
    # j, s ..
end

@enum WaveEquationFormat schrodinger # ALSO dirac and scarel for scalar_relativistic

"""Compute eigenvalue of the hydrogen-like atoms for bound of solution"""
function hydrogen_compute(Z::Int, qn::QuantumNumber, wqf::WaveEquationFormat)
    if wqf === schrodinger
        ev = -Z ^ 2 / (2 * qn.n ^ 2)
    end

    ev
end
    

end
