<img src="https://github.com/unkcpz/PseudopotentialGenerator.jl/blob/main/misc/logo/_light.png?raw=true" alt="pgen logo" height="100px" />

# PseudopotentialGenerator.jl

| **Build Status**                                |  **License**                     |
|:----------------------------------------------- |:-------------------------------- |
| [![][ci-img]][ci-url] [![][ccov-img]][ccov-url] | [![][license-img]][license-url]  |

[ci-img]: https://github.com/unkcpz/PseudopotentialGenerator.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/unkcpz/PseudopotentialGenerator.jl/actions

[ccov-img]: https://codecov.io/gh/unkcpz/PseudopotentialGenerator.jl/branch/main/graph/badge.svg?token=2KH3oPQm9E
[ccov-url]: https://codecov.io/gh/unkcpz/PseudopotentialGenerator.jl

[license-img]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-url]: https://github.com/unkcpz/PseudopotentialGenerator.jl/blob/main/LICENSE

`PseudopotentialGenerator.jl` generates pseudopotential for plane-wave DFT.

## Usage

Generating a pseudopotential requires first having reference states (from real all-electrons potential) and then construct the pseudopotential that can reproduce same scattering behavior.

The [examples](https://github.com/unkcpz/PseudopotentialGenerator.jl/tree/main/examples) folder includes scripts of solving atomic Schrödinger equation and pseudizing to get pseudopotentials.

To solve the atomic Schrödinger equation under DFT framework, simply define the ionic configuration and call `self_consistent_field` function. 

```julia
Z = 6
mesh = Mesh(1e-7, 50.0, 2.7e+6, 10000);
orbs = [Orbital(Z, 1, 0, 2), Orbital(Z, 2, 0, 2), Orbital(Z, 2, 1, 2)];

res = self_consistent_field(
    Z,
    mesh,
    orbs,
    mixing_beta = 0.2,
    abstol = 1e-7,
    maxiters_scf = 200,
    perturb = true,
)

ae_info = (ε_lst = res.ε_lst, ϕs = res.ϕs, vae = res.v_tot, occs = res.occs, xc = res.xc);
```

The `ae_info` contains eigenstates and self consistent potential.

The pseudopotential is then generated to have pseudo-states that meet the restriction with respect to the all-electrons reference.
The [Troullier-Martins (TM)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.43.8861) is a one of the popular scheme to construct pseudopotential.
To generate a TM type pseudopotential, call `pseudize` function,

```julia
v_pspot, ϕ_ps = pseudize(ae_info, mesh, rc; method = :TM, kbform = false)
```

The function returns the psudopotential and all pseudo-wavefunctions.

It is essential to check whether the generated pseudopotential is good enough for solid-state DFT calculations.
One of the prevalent metric is logarithmic derivative (usually `arctan` of the derivative) of wavefunctions which reveal the discripancies between the all-electrons wavefunctions and pseudo-wavefunctions in both eigenstates and energies of scattering states.
The arctan logarithmic derivative is computed for every angular momentum by calling `compute_atanld` in the given potential,

```julia
# For all-electrons potential
compute_atanld(l, Z, mesh, ae_info.vae, rcx, window = [-10.0, 10.0], δ = 0.05)

# For pseudo-potential which is semi-local (dependent on angular momentum l)
compute_atanld(l, Z, mesh, v_pspot.v[nl], rcx, window = [-10.0, 10.0], δ = 0.05)
```

## Basic theory and package construction

You will find presentian that covers the thoretical background and the construction detials of this package.

[DFT in atom, pseudopotential generation and its Julia implementation](https://raw.githubusercontent.com/unkcpz/PseudopotentialGenerator.jl/refs/heads/main/misc/PGEN_atomic_DFT_in_Julia.pdf)

## For developers

To run tests, functions from [DFTATOM](https://github.com/certik/dftatom) are wrappered as reference which need to compiled first.

```bash
make -C deps all
```

Then run tests

```bash
julia --project=@. -e 'using Pkg; Pkg.test()'
```

## TODOs

- [x] Solving atomic DFT in radial coordination
- [x] Migrate from self-crafted ODE solver to purely using `NonlinearSolve.jl` and `OrdinaryDiffEq.jl` as solver.
- [x] pseudize using TM.
- [ ] pseudize using RRKJ.
- [ ] pseudize using BHS.
- [x] compute logarithmic derivative plots.
- [ ] Kleinman-Bylander form and its logarithmic derivative by solving integ-ODE.
- [ ] full relativistic PP
- [ ] scalar relativistic PP
- [ ] Orbitals defined and constructed from atomic levels representation.
- [ ] Logger for SCF and eigensolver
- [ ] Support for semi-core states
- [ ] Support for ONCV NC pseudopotentials
- [ ] Support for US/PAW pseudopotentials
- [ ] meta-GGA PP generation
- [ ] ? Scattering states with good match
