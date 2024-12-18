module PseudopotentialGenerator

# Reference implementation for testing
include("DFTATOM.jl")
export FPGEN    # FPGEN is the FORTRAN wrapper
export DFTATOM_JL   # DFTATOM_JL is the Julia implementation

# Julia implementation
include("common.jl")
export SPEED_OF_LIGHT, Orbital

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
export poisson_outward

include("pseudize.jl")
export pseudize

include("ld.jl")
export compute_ld, compute_atanld

macro ignore(args...) end

# For LSP working at example scripts
@ignore include("../examples/logder.jl")
@ignore include("../examples/continuity_envolve.jl")
@ignore include("../examples/numerical_stability.jl")

# at tests
@ignore include("../test/runtests.jl")

end

