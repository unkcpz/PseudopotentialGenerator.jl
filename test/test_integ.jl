@testset "ode1d Newton-Cotes" begin
    x = range(0, 3, length=20)
    y = sin.(x)
    xp = cos.(x)
    for method in [:trapz1, :trapz3, :trapz5, :trapz7, :simpson, :adams]
        f_method = Symbol("f_", method)
        @test integrate(y, xp, method=method) ≈ FPGEN.integrate(y, xp, method=f_method) atol = 1e-12
    end
end

"""
    analytical_solution_damped(x0, v0, zeta, omega, t)

Define the analytical solution for a damped harmonic oscillator with initial conditions x0 and v0.

Using damped harmonic oscillator to test the integration of 2nd-ored ode.
The second order ODE is:
y'' + 2 * ζ * ω * y' + ω² * y = 0

convert to 2 1st-order ODEs:
x' = v
v' = - ω² * x - 2 * ζ * ω * v

where x = y and v = y'

The analytical solution is:
x(t) = exp(-ζ * ω * t) * (x0 * cos(ω_d * t) + ((v0 + ζ * ω * x0) / ω_d) * sin(ω_d * t)
v(t) = -ζ * ω * exp(-ζ * ω * t) * (x0 * cos(ω_d * t) + ((v0 + ζ * ω * x0) / ω_d) * sin(ω_d * t) +
       exp(-ζ * ω * t) * (-x0 * ω_d * sin(ω_d * t) + (v0 + ζ * ω * x0) * cos(ω_d * t))

where ω_d = ω * √(1 - ζ²)

Note: The problem is not the same as the radial Schrödinger equation as we have. The difference
is that the in one of two ODEs, the coefficient of the second derivative is depend on r (or t in the case of the damped harmonic oscillator).
"""
function analytical_solution_damped(x0, v0, zeta, omega, t)
    omega_d = omega * sqrt(1 - zeta^2)
    x_t = exp(-zeta * omega * t) * (x0 * cos(omega_d * t) + ((v0 + zeta * omega * x0) / omega_d) * sin(omega_d * t))

    e_term = exp(-zeta * omega * t)
    cos_term = cos(omega_d * t)
    sin_term = sin(omega_d * t)

    v_t = -zeta * omega * e_term * (x0 * cos_term + ((v0 + zeta * omega * x0) / omega_d) * sin_term) +
          e_term * (-x0 * omega_d * sin_term + (v0 + zeta * omega * x0) * cos_term)

    return x_t, v_t
end

@testset "ode1d rk4" begin
    mesh = Mesh(0.0, 0.5, 1.0, 100)
    r = mesh.r
    y0 = [1.0, 0.0]
    ζ = 0.1
    ω = 2.0 * π
    C1 = C1mid = fill(-ω^2, length(mesh.r))
    C2 = C2mid = fill(-2 * ζ * ω, length(mesh.r))
    max_val = 10.0    # Upper bound value of the solution

    y1_dftatom, y2_dftatom, imax = FPGEN.rk4_integrate(mesh.r, y0, C1, C2, C1mid, C2mid, max_val)

    @test imax == length(mesh.r)

    # Expected analytical solution
    for idx in [1, 50, 100]
        y1_expected, y2_expected = analytical_solution_damped(y0[1], y0[2], ζ, ω, r[idx])

        @test y1_dftatom[idx] ≈ y1_expected atol = 1e-5
        @test y2_dftatom[idx] ≈ y2_expected atol = 1e-5
    end

    # Test my implementation of rk4 give the same results as dftatom (direct translate)
    y0 = [1.0, 0.0]
    y1, y2, imax = rk4_integrate(mesh.r, y0, C1, C2, C1mid, C2mid, max_val)

    @test imax == length(mesh.r)

    for idx in [1, 50, 100]
        @test y1[idx] ≈ y1_dftatom[idx] atol = 1e-12
        @test y2[idx] ≈ y2_dftatom[idx] atol = 1e-12
    end

    # Test the case where C1 and C2 are constant
    # Test using the function f1 and f2

    # f1 = dy1/dr -> y2
    # f2 = dy2/dr -> C1 * y1 + C2 * y2
    f1(y1, y2, r) = y2
    f2(y1, y2, r) = -ω ^ 2 * y1 - 2 * ζ * ω * y2
    y0 = [1.0, 0.0]
    y1, y2, imax = rk4_integrate(f1, f2, mesh.r, y0, max_val)

    @test imax == length(mesh.r)

    for idx in [1, 50, 100]
        @test y1[idx] ≈ y1_dftatom[idx] atol = 1e-12
        @test y2[idx] ≈ y2_dftatom[idx] atol = 1e-12
    end
end

@testset "ode1d rk4 with r depended C1/C2" begin
    # Test the C1/C2 can be a function of r
    mesh = Mesh(0.0, 0.5, 1.0, 100)
    max_val = 10.0    # Upper bound value of the solution
    r = mesh.r
    ζ = 0.1
    ω = 2.0 * π
    y0 = [1.0, 0.0]
    C1 = -ω^2 * r
    C2 = -2 * ζ * ω * r
    C1mid = midpoints(C1, r)
    C2mid = midpoints(C2, r)

    y1_dftatom, y2_dftatom, imax = FPGEN.rk4_integrate(mesh.r, y0, C1, C2, C1mid, C2mid, max_val)

    f1(y1, y2, r) = y2
    f2(y1, y2, r) = -ω ^ 2 * r * y1 - 2 * ζ * ω * r * y2
    y0 = [1.0, 0.0]
    y1, y2, imax = rk4_integrate(f1, f2, mesh.r, y0, max_val)

    @test imax == length(mesh.r)

    for idx in [1, 50, 100]
        @test y1[idx] ≈ y1_dftatom[idx] atol = 1e-12
        @test y2[idx] ≈ y2_dftatom[idx] atol = 1e-12
    end

end

@testset "rschored" begin
    V(r) = 1 / r # Coulomb like potential

    for l in [0, 1]
        for Z in [1, 1]
            for E in [-0.5]
                mesh = Mesh(1e-6, 50.0, 1e+9, 500)

                Q_expect, P_expect, imax_expect = FPGEN.schroed_outward_adams(l, Z, E, V, mesh.r, mesh.rp)
                Q, P, imax = schroed_outward_adams(l, Z, E, V, mesh.r, mesh.rp)
                @test imax == imax_expect
                for idx in [1, 50, 100]
                    @test Q[idx] ≈ Q_expect[idx] atol = 1e-12
                    @test P[idx] ≈ P_expect[idx] atol = 1e-12
                end

                Q_expect, P_expect, imin_expect = FPGEN.schroed_inward_adams(l, E, V, mesh.r, mesh.rp)
                Q, P, imin = schroed_inward_adams(l, E, V, mesh.r, mesh.rp)
                @test imin == imin_expect
                for idx in 400:2:405
                    @test Q[idx] ≈ Q_expect[idx] atol = 1e-12
                    @test P[idx] ≈ P_expect[idx] atol = 1e-12
                end

            end
        end
    end
end



@testset "poisson" begin
    # Test the Poisson solver
    mesh = Mesh(1e-7, 10.0, 1.0, 100)

    # wave function in infinite square well potential
    a = 5
    n = 2
    ψ = sqrt(1/a) * sin.(n * π * mesh.r / a)

    # charge density
    ρ = abs2.(ψ)

    # Poisson solver
    vh_expected = FPGEN.rpoisson_outward_pc(ρ, mesh.r, mesh.rp)
    vh = rpoisson_outward_pc(ρ, mesh)
    @test vh[1:10] ≈ vh_expected[1:10] atol = 1e-12
end