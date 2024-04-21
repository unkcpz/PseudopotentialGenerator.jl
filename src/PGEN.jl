module PGEN

# FORTRAN wrapper
include("FPGEN.jl")
export FPGEN

# Julia implementation
include("common.jl")
export SPPED_OF_LIGHT, Orbital

include("ode1d.jl")
export integrate, rk4_integrate

include("mesh.jl")
export Mesh, mesh_exp_deriv2, midpoints, dfdr, d2fdr2

include("rschroed.jl")
export schroed_outward_adams, schroed_inward_adams, sch_outward, sch_inward

include("reigen.jl")
export solve_radial_eigenproblem, find_ctp

include("potential.jl")
export coulomb_potential, thomas_fermi_potential
export SemiLocalPotential, KBFormPotential, Potential

include("scf.jl")
export self_consistent_field

include("rpoisson.jl")
export rpoisson_outward_pc, poisson_outward

include("pseudolize.jl")
export pseudolize

include("ld.jl")
export compute_ld, compute_atanld

end
