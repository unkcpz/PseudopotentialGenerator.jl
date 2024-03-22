using Libdl

libDFTATOM = joinpath(@__DIR__, "../../deps/libdftatom.so")

# ODE 1D solver
include("ode1d.jl")
export integrate