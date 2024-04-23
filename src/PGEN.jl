module PGEN

# Reference implementation for testing
include("DFTATOM.jl")
export FPGEN    # FPGEN is the FORTRAN wrapper
export DFTATOM_JL   # DFTATOM_JL is the Julia implementation

# Julia implementation
include("common.jl")
export SPPED_OF_LIGHT, Orbital

include("integ1d.jl")
export integrate

include("mesh.jl")
export Mesh, mesh_exp_deriv2, midpoints, dfdr, d2fdr2

include("ode_sch.jl")
export sch_outward, sch_inward

include("eigensolver.jl")
export solve_radial_eigenproblem, find_ctp

include("potential.jl")
export coulomb_potential, thomas_fermi_potential
export SemiLocalPotential, KBFormPotential, Potential

include("scf.jl")
export self_consistent_field

include("ode_poisson.jl")
export rpoisson_outward_pc, poisson_outward

include("pseudolize.jl")
export pseudolize

include("ld.jl")
export compute_ld, compute_atanld

end
