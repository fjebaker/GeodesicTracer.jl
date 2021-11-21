# Second order integration

```@meta
CurrentModule = GeodesicTracer
```

The philosophy behind the second order methods is to permit the computation of geodesics in generic spacetimes, via the geodesic equation

```math
\frac{\text{d}^2 x^\mu}{\text{d} \lambda^2} 
    + \Gamma^{\mu}_{\phantom{\mu}\nu\sigma}
    \frac{\text{d}x^\nu}{\text{d} \lambda}
    \frac{\text{d}x^\sigma}{\text{d} \lambda}
= 0,
```

where $x^\mu$ is a position four-vector, $\Gamma^{\mu}_{\phantom{\mu}\nu\sigma}$ are the Christoffel symbols of the second kind, and $\lambda$ the affine parameter describing the curve.

The above can be solved as a second order ODE, subject to an initial position

```math
x^\mu = \left(t_0, x_0, y_0, z_0 \right),
```

where all components are known. Note, any coordinate system can be used in place of $x_0$, $y_0$. $z_0$, but at the moment the code is configured towards using spherical coordinates -- this will be relaxed as additional coordinate systems are used.

The integration also requires an initial velocity 

```math
\frac{\text{d}x^\mu}{\text{d} \lambda} = \left( \dot{t}_0, \dot{x}_0, \dot{y}_0, \dot{z}_0 \right),
```

where the dot refers to the derivative with respect to $\lambda$. In general, $\dot{x}_0, \dot{y}_0, \dot{z}_0$ are known, and $\dot{t}_0$ is computed to constrain the geodesic's behaviour. That is, $\dot{t}_0$ is chosen such that the relation

```math
g_{\mu\nu} \frac{\text{d}x^\mu}{\text{d} \lambda} \frac{\text{d}x^\nu}{\text{d} \lambda} = \kappa^2,
```

with the metric tensor $g_{\mu\nu}$, and where $\kappa$ can be related to mass. For null geodesics, $\kappa^2 = 0$ always, used to calculate the trajectory of light in a given spacetime. Time-like trajectories depend on $\text{sign}(g_{\mu\nu})$.

The minimal configuration for a new spacetime is to have the geodesic equation and a constraint equation defined. [DiffGeoSymbolics.jl](https://github.com/astro-group-bristol/DiffGeoSymbolics.jl) will provide Julia tools for doing this automatically, and will be at some point included in this library. [ComputedGeodesicEquations.jl](https://github.com/astro-group-bristol/ComputedGeodesicEquations.jl) provides SageMath computed symbolic expression, which have been transformed into (relatively) optimised Julia functions for different spacetimes.

## Defining integration problems

The dispatch for different spacetimes is controlled by the parametric parameter of:

```@docs
GeodesicParams
```

The two (three, if counting the in-place variant) equations that must be defined are then [`secondorder_rayintegrator`](@docs), and [`null_constrain`](@docs) (*todo* allow non-null geodesics).

## Integration methods

These are currently defined in the default library. Other libraries may in the future extend this for specific problems.

```@docs
secondorder_rayintegrator
```

## Solution handling

The solutions return contain vectors where the first 4 elements are the velocity components, and the last 4 are the position. In practice, some solution `sol` can be evaluated at time `λ` and decomposed via

```julia
u = sol(λ)
v, x = @view(u[1:4]), @view(u[5:8])
```

or, via the internal representation:

```julia
u = sol(λ)
v, x = u.x
```