#using Revise
using PseudopotentialGenerator
using Test

include("test_legacy.jl")
include("test_mesh.jl")
include("test_scf.jl")
include("test_integ.jl")
include("test_eigensolver.jl")
include("test_pseudize.jl")
