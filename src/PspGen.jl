module PspGen

using Interpolations
using DifferentialEquations
# using Printf

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

"""Compute eigenvalue of the hydrogen-like atoms for bound of solution
Z for atomic number, n for principal quantum number"""
function compute_hydrogen(Z::Int, n::Int)
    ev = -Z ^ 2.0 / (2 * n ^ 2)

    ev
end

# compute min of eigenvalue from potential (XXX)
function compute_emin(Z::Int, n::Int)
    0.98 * compute_hydrogen(Z, n)   # a little bit small than hydrogen 1s.
end

# compute max of eigenvalue from potential
function compute_emax()
end

# compute range of every eigenvalue
function eigensolver_bracket()
end

struct Potential
    xs::Vector{Float64}
    ys::Vector{Float64}

    Potential(xs, ys) = length(xs) != length(ys) ? error("dims not matched") : new(xs, ys)
end

# boundary condition of point close to origin r0
# return tuple for u, u_p
function bc_origin(r0::Float64, Z::Int, l::Int, ε::Float64, p::Potential)
    # SCHRODINGER
    # For a potential v(r) = vp(r) - Z/r
    # -Z / r is the external part and vp(r) is the local (v_hatree + v_xc) potential at r.
    # The general solutions of the equations when r->0 are of the form:
    #
    #   u(r) = r^(s) (a0 + a1 r + a2 r^2 + ...)
    #   u_p(r) = r^(s-1) (b0 + b1 r + b2 r^2 + ...)
    # 
    # Schrodinger equation:
    #   s = l
    #   a0 = 0
    #   a1 = 1 (set to 1, only scale the function not influence eigenvalue)
    #   a2 = - Z * a1 / (l + 1)
    #   a3 = [-Z * a2 + (vp - e) a1] / (2l + 3)
    #   b0 = l a1
    #   b1 = (l + 1) a2
    #   b2 = (l + 2) a3
    vp0 = v(r0, p) + Z / r0
    s = l
    a0 = 0
    a1 = 1
    a2 = -Z * a1 / (s + 1)
    a3 = (-Z * a2 + (vp0 - ε)) / (2 * s + 3)
    b0 = s
    b1 = (s + 1) * a2
    b2 = (s + 2) * a3

    u = r0 ^ (s-1) * (a0 + a1 * r0 + a2 * r0 ^ 2 + a3 * r0 ^ 3)
    u_p = r0 ^ (s-1) * (b0 + b1 * r0 + b2 * r0 ^ 2)

    u, u_p
end

# return tuple for u, u_p at very large r, V(r) -> 0 and it is the a ode solved analytically.
function bc_infinity(ε::Float64)
    FUNC_INF = eps(Float64) # function value at r=inf

    # r->∞: u = exp(-√(-2ε))
    u = FUNC_INF
    u_p = -√(-2ε) * u

    u, u_p
end

# get potential at r using irregular interplate
function v(r::Float64, p::Potential)
    interp = interpolate(p.xs, p.ys, FritschCarlsonMonotonicInterpolation())
    
    interp(r)
end

function SE_ODE!(du, u, p, r)
    # SCHRODINGER 2nd order ODE
    #
    #                    d R
    # y_1 = R     y_2 =  ---
    #                    d r
    #
    # f_1 = y_2
    #
    #         2        l(l+1)
    # f_2 = - - y_2 + [------ + 2(v_ks - e)] y_1
    #         r         r^2       
    # l, ε, pot = p
    # rinv = 1 / r
    # r2i = rinv * rinv
    # pot -> Potential
    v(r) = -1/r
    f(r, e) = v(r) - e
    du[1] = u[2]
    du[2] = f(r, p) * u[1]
end

"""
The log der diff at turning point rm at given eigenvalue e.
"""
function ldd_rm(ε::Float64)
    rm = -1/ε   # f(r, e) = -1/r - e

    prob(u0, rspan, e) = ODEProblem(SE_ODE!, u0, rspan, e)

    ODE_INF = 1.0e-12   # The practical value of function at INF.
    
    rmin = 1.0e-8  # the smallest grid x value of potential
    a1 = 1 
    a2 = -1/2
    a3 = 1/12 - ε/6

    r = rmin
    a1 = 1
    a2 = -2
    a3 = 1/3
    a4 = -1/36
    bc = a1 * r + a2 * r^2 + a3 * r^3
    bcp = a1 + 2 * a2 * r + 3 * a3 * r^2 + 4 * a4 * r^3
    bc_zero = [bc, bcp]

    # emax = -1.0e-12 # either the smallest possible eigenvalue or a very small value
    rmax = -log(ODE_INF) / sqrt(-1*ε)   # The solution when set V as 0. use the smallest ε if there are more.

    bc_inf = [ODE_INF, -√(-1/ε) * ODE_INF]

    # outward from origin, and inward from inf
    ldd = logder_diff(prob, ε, bc_zero, bc_inf, rmin, rmax, rm)
    ldd
end

function logder_diff(
    prob, 
    ε::Float64, 
    bc_zero, 
    bc_infinity, 
    rmin::Float64, 
    rmax::Float64, 
    rm::Float64;
    algo=Tsit5(), 
    abstol::Float64=1e-12,
)
    outward = solve(prob(bc_zero, (rmin, rm), ε), algo, abstol=abstol, save_everystep=false)
    inward = solve(prob(bc_infinity, (rmax, rm), ε), algo, abstol=abstol, save_everystep=false)

    ldd = outward.u[end][2] / outward.u[end][1] - inward.u[end][2] / inward.u[end][1]
    ldd
end

end