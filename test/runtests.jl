using Revise
using Test
using PGEN

include("test_legacy.jl")
include("test_mesh.jl")
include("test_scf.jl")
include("test_integ.jl")
include("test_eigensolver.jl")
include("test_pseudolize.jl")


mesh = Mesh(1e-6, 50.0, 1e+7, 5000)

Z = 10
l = 1

P, Q, imax = DFTATOM_JL.schroed_outward_adams(l, Z, E, V, mesh.r, mesh.rp)
P_, Q_, imax_ = sch_outward(l, Z, E, V, mesh.r, mesh.rp)

using Plots

plot(mesh.r, P, label="P")
plot!(mesh.r, P_, label="P_")