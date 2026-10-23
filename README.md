# NeutrinoMassModels.jl

Models of neutrino masses for [GLoBES.jl](https://github.com/arnavk23/GLoBES.jl): particle models written as Mathematica-style model files, their tree-level spectra, and neutrino mass mechanisms that turn model parameters into oscillation parameters. The package was split out of GLoBES.jl, which keeps the oscillation physics.

```julia
using Pkg
Pkg.develop(path = "/path/to/NeutrinoMassModels.jl")   # not registered yet
using NeutrinoMassModels, GLoBES
```

## Particle models

`load_particle_model` reads a model directory of Mathematica-style files (the model file, `particles.m`, `parameters.m`, `SPheno.m`); `BSMModel` and `add_gauge!`, `add_fermion!`, `add_scalar!`, `add_lagrangian!` build one in Julia.

```julia
scoto = load_particle_model("Scotogenic/")
check_model(scoto)                # gauge and global invariance, mass dimension of each term
s = spectrum(scoto, values)       # tree-level scalar and Majorana masses, tadpoles solved
write_leshouches("LesHouches.in", scoto, values)
```

## Neutrino masses

`ScotogenicLoop` (E. Ma, hep-ph/0601225) and `TypeISeesaw` give the light neutrino mass matrix of a model, and `casas_ibarra` the Yukawa couplings that reproduce oscillation data:

```julia
mech = ScotogenicLoop()           # parameter names as in the scotogenic model file
values[:Yn] = casas_ibarra(mech, values, nufit_params("6.1"); lightest = 1e-3)
neutrino_spectrum(mech, values)   # masses, angles, phases, mββ, Σm
br_lfv(mech, values, :μ, :e)      # BR(μ → eγ); MEG II: 1.5e-13
```

`latex(model)`, `latex(check_model(model))`, `latex(spectrum(...))` and `latex(neutrino_spectrum(...))` write the results for a GLoBES.jl `TexReport`.

Tests: `julia --project=. test/runtests.jl`.
