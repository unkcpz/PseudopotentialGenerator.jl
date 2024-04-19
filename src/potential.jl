abstract type Potential end

struct SemiLocalPotential <: Potential
    v::Dict{NamedTuple{(:n, :l), Tuple{Int64, Int64}}, Vector{Float64}}
end

struct KBFormPotential <: Potential
    v_local::Vector{Float64}
    ekb
    proj_kb
end

# TODO: ? also make it a Potential type?
function coulomb_potential(Z::Int64, r::Float64)::Float64
    return -Z / r
end

function thomas_fermi_potential(Z::Int64, r::Float64)::Float64
    # from DFTATOM
    # Z_eff(x) = Z * phi(x), where phi(x) satisfies the Thomas-Fermi equation:
    #   phi'' = phi**(3/2) / sqrt(x)
    # with boundary conditions:
    #   phi(0)  = 1
    #   phi(oo) = 0
    # There is no analytic solution, but one can solve this approximately. We use:
    # http://arxiv.org/abs/physics/0511017

    x = r * (128 * Z / (9*π^2))^(1/3)

    α = 0.7280642371
    β = -0.5430794693
    γ = 0.3612163121

    Z_eff = Z * (1 + α * sqrt(x) + β * x * exp(-γ * sqrt(x)))^2 * exp(-2 * α * sqrt(x))

    # always cut the potential to be greater than 1
    # so the eigenvalues are negative
    Z_eff = max(1, Z_eff)
    V = -Z_eff / r
    V
end