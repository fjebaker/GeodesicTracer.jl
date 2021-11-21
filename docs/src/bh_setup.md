# Setup & Configuration

```@meta
CurrentModule = GeodesicTracer
```

## Spacetime setup

All geodesic problems are configured via the Black Hole Setup (`BHSetup`) structure:

```@docs
BHSetup
rescale!
```

## Integration setup

The main function used to calculate geodesics is:

```@docs
calcgeodesic
```

In the case where specific control of an integration problem is needed, GeodesicTracer.jl presents configuration classes for serial and parallel computation:

```@docs
IntegratorConfig
SingleParams
ParallelParams
```

For `ParallelParams`, a problem function must be defined that updates the initial conditions for each ensemble: 

```@docs
makeprobfunc
```


These can be instances and used with:

```@docs
integrategeodesic
```

## Integration chart

All geodesics are light-like (*todo*: time-like), starting by default at $λ=0$, $r=1000$ orientated towards the black hole.

The integrator defines a chart 

```@docs
chartbounds
```

which interrupts the integration when the geodesic leaves the local spacetime.

### Configuration

We distinguish between two common problem types

- those where the geodesic paths are desired
- those where the collision points of the geodesics are desired

This is done for performance; the memory overhead of the later can be far reduced by not storing intermediate values in the integration process.

Internally, these are differentiated through the keyword parameter `save_geodesics` in [`calcgeodesic`](@ref) or passed to an [`IntegratorConfig`](@ref):
