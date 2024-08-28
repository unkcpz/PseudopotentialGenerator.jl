# Speed of light in atomic units
const SPEED_OF_LIGHT = 137.035999037

# Hydrogen-like atom energy: E = - Z^2 / (2.0 * n^2)
function hydrogen_like_energy(n::Int64, Z::Int64)::Float64
    -Z^2 / (2.0 * n^2)
end

struct Orbital
    Z::Int64
    n::Int64
    l::Int64
    f::Int64
    emin::Float64
    emax::Float64

    function Orbital(Z::Int64, n::Int64, l::Int64, f::Int64)
        emin = hydrogen_like_energy(n, Z) * 1.1 # 10% for robustness
        emax = 10.0
        new(Z, n, l, f, emin, emax)
    end
end
