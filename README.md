# GeodesicTracer.jl

## Setup

Until the module is registered with the Julia package manager, the easiest way to use the project is just to clone the source code.

To setup from a fresh Julia install:

- clone repository
```bash
https://github.com/astro-group-bristol/GeodesicTracer.jl && cd GeodesicTracer.jl
```

- initialize the project in the Julia REPL with desired number of threads (e.g. 8)
```bash
julia --project=. --threads=8
```

- install dependencies and pre-compile
```julia
julia> ]initialise
```

- (optional) run tests
```julia
julia> ]test
```

## Example use

```julia
using GeodesicTracer
using Plots

# setup
s = BHSetup(a=1.0)

# set impact parameter ranges
α_range = (0.1, 4.0)
β = 0.1

# number of geodesics to calculate in that range
num = 100


simsol = @time calcgeodesic(α_range, num, β, s)

# plot ISCO
test_plot = plot(
    _ -> GeodesicTracer.rms(s), 
    0:0.01:2π, 
    line = (:dot, 2),
    proj = :polar,
    color = :black,
    legend = false
)

# plot event horizon
plot!(
    _ -> R₀(s), 
    0:0.01:2π, 
    lw = 4, 
    color = :black
)

# plot geodesics
plot!(
    simsol,
    vars=(4, 2) # choose phi and r as components
)

# show figure up to r = 30 MG
plot(test_plot, range=(0, 30))
```

To produce a geometric disk render

```julia
using GeodesicTracer
using Plots

s = BHSetup(a=1.0)

d = GeometricDisk(
    r_inner=GeodesicTracer.rms(s), 
    inclination=deg2rad(80.0)
)

image = @time render(d, s)

ht = heatmap(image, clims=(0.0, 2.0))
ht
```

For the shadow of a black hole, we use a custom value function which returns just the coordinate time
```julia
using GeodesicTracer
using Plots

s = BHSetup(a=1.0)

const time = ValueFunction(
    (val, λ, u, p, d) -> λ
)

image = @time render(s, valuefunc=time)

ht = heatmap(img, clim=(497.0, 510.0))
ht
```


## Documentation

To build the documentation
```bash
julia --project=. docs/make.jl
```
which will create and populate the directory `docs/build`, containing a static webpage with auto documentation.

Note, this can take a while, since `Documenter.jl` needs to load full Julia environment to properly inspect the code.