# Configuration

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

In the case where specific control of an integration problem is needed, GeodesicTracer.jl presents configuration classes for serial and parallel computation:

```@docs
IntegratorConfig
SingleParams
ParallelParams
```

These can be instances and used with:

```@docs
integrategeodesic
```