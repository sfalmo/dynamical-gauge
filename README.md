# Dynamical gauge invariance of statistical mechanics

This is a proof-of-concept implementation of a 1D nonequilibrium molecular dynamics simulation which computes gauge correlation functions via automatic and finite initial-state time differentiation as described in:

**Dynamical gauge invariance of statistical mechanics**  
*Johanna Müller, Florian Sammüller, and Matthias Schmidt, to be published.*

## Instructions

A recent version of [Julia](https://julialang.org/downloads/) needs to be installed on your system.
Launch the Julia interpreter within this directory and type `]` to enter the package manager.
Activate the environment and install the required packages as follows:

```julia
activate .
instantiate
```

Type backspace to exit the package manager.

Fetch pregenerated data (~2GB) and reproduce the plots from the manuscript:

```julia
using Downloads, Tar
Downloads.download("https://www.staff.uni-bayreuth.de/~bt306964/dynamical-gauge/data.tar", "data.tar")
Tar.extract("data.tar")
include("plot_paper.jl")
```

Run a custom nonequilibrium simulation:

```julia
include("main.jl")
```

Inspect the sampled gauge correlation functions:

```julia
include("plot_all_histograms.jl")
```
