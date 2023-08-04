
"""
    lschfb solve the radial eigenvalue problem for a given potential.
        wrapping and calling lschfb from oncvpsp.
"""
LIBONCV_PATH = joinpath(@__DIR__, "../../deps/liboncv.so")

# Thomas-Fermi potential
export tf
include("tfapot.jl")

# Hartree potential
#include("hartree.jl")
#include("lsch.jl")
