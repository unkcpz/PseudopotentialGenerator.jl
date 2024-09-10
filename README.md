<img src="https://github.com/unkcpz/PGEN.jl/blob/main/misc/logo/PGEN_light.png?raw=true" alt="pgen logo" height="100px" />

# PGEN.jl

| **Build Status**                                |  **License**                     |
|:----------------------------------------------- |:-------------------------------- |
| [![][ci-img]][ci-url] [![][ccov-img]][ccov-url] | [![][license-img]][license-url]  |

[ci-img]: https://github.com/unkcpz/PGEN.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/unkcpz/PGEN.jl/actions

[ccov-img]: https://codecov.io/gh/unkcpz/PGEN.jl/branch/main/graph/badge.svg?token=2KH3oPQm9E
[ccov-url]: https://codecov.io/gh/unkcpz/PGEN.jl

[license-img]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-url]: https://github.com/unkcpz/PGEN.jl/blob/main/LICENSE

PGEN.jl generates pseudopotential for plane-wave DFT.

## Usage

The [examples](https://github.com/unkcpz/PGEN.jl/tree/main/examples) folder provide scripts include solving atomic Schrödinger equation and pseudolize to get pseudopotentials.

## For developers

To run tests, some functions from [DFTATOM](https://github.com/certik/dftatom) is wrappered as reference which need to compiled first.

```bash
make -C deps all
```

Then can run tests by 

```bash
julia --project=@. -e 'using Pkg; Pkg.test()'
```

## TODO

- [x] Solving atomic DFT in radial coordination
- [x] Migrate from self-crafted ODE solver to purely using `NonlinearSolve.jl` and `OrdinaryDiffEq.jl` as solver.
- [x] pseudolize using TM.
- [ ] pseudolize using RRKJ.
- [ ] pseudolize using BHS.
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
