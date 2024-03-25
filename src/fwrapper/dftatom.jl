using Libdl

libDFTATOM = joinpath(@__DIR__, "../../deps/libdftatom.so")

# ODE 1D solver
include("ode1d.jl")
export integrate, rk4_integrate

# Mesh (exponential grid)
include("mesh.jl")
export mesh_exp, mesh_exp_deriv, mesh_exp_deriv2

# radial Schrödinger equation outward and inward integration
include("rshroed.jl")