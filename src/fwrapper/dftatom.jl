using Libdl

libDFTATOM = joinpath(@__DIR__, "../../deps/libdftatom.so")

# ODE 1D solver
include("ode1d.jl")
export integrate

# Mesh (exponential grid)
include("mesh.jl")
export mesh

# radial Schrödinger equation outward and inward integration
include("rshroed.jl")