# First order integration

```@meta
CurrentModule = GeodesicTracer
```

First order geodesic methods rely on being able to solve an ODE system separably, which requires linearising the geodesic equation such as via the Hamilton-Jacobi method. This is often non-trivial, and introduces new manifolds and embeddings that require the integration method to do something extra, such as keeping track of additional changing quantities.

As a consequence, the first order methods require special treatment from the integrator perspective, and tailored methods specific to a particular solution.

Below are the currently implemented solutions in this library (which is currently only one).

## Carter's linearized solutions in Boyer-Lindquist coordinates

Carter (1968) derived a fourth constant of motion from the Boyer-Lindquist coordinates describing the Kerr spacetime, which allowed the 2nd order geodesic equations to be linearized[^1].

[^1]: B. Carter (1968), *Global Structure of the Kerr Family of Gravitational Fields*, Phys Rev, **174**, 5.

The Kerr metric in the Boyer-Lindquist coordinates is written

```math
\begin{aligned}
\text{d}s^2 = 
& - \left( 1 - \frac{2Mr}{\Sigma} \right) \text{d}t^2 \\
& - \left( \frac{4 Mar \sin^2}{\Sigma} \right) \text{d}t \text{d}\phi \\
& + \left( \frac{\Sigma}{\Delta} \right) \text{d}r^2 \\
& + \Sigma \text{d}\theta^2 \\
& + \left( 
        r^2 + a^2 + \frac{2 M a^2 r \sin^2 \theta}{\Sigma}
    \right) \sin^2 \theta \text{d}\phi^2,
\end{aligned}
```

with:

```@docs
Σ
Δ
```

The metric structure allowing these to be changed is

```@docs
CarterBoyerLindquist
```

## Equations of motion

The following separable ODE system describes a photon trajectory in the Kerr spacetime:

```@docs
Σδt_δλ
Σδr_δλ
Σδθ_δλ
Σδϕ_δλ
```

with potentials

```@docs
Vr
Vθ
```

In terms of the solver we implement, these four equations are wrapped in
```@docs
δ 
δ!
```

The dispatch calling these:

```@docs
rayintegrator(u, p::CarterGeodesicParams, λ)
rayintegrator!(du, u, p::CarterGeodesicParams, λ)
```

## Constants of motion

The constants of motion are:

```@docs
T
S
LQ
```

In `LQ`, we additionally have:

```@docs
A
```


## Properties of the spacetime

The event horizon and marginally stable / innermost stable circular orbit (ISCO) is given by: 

```@docs
R₀
rms
Z₁
Z₂
```

## Impact parameter mapping

The impact parameters $\alpha$ and $\beta$ are mapped into observer angles via:

```@docs
sinΦsinΨ
```


## Implementation details

The sign of the potentials [`Vr`](@ref), [`Vθ`](@ref) is tracked via `DiscreteCallback` functions. This allows the potentials to be checked multiple times per time step (depending on the solver), but updated only once per time step.

This is achieved via

```@docs
is_radial_pot_negative
is_angular_pot_negative
flip_radial_sign
flip_angular_sign
```

Ensemble renderings create new integration problems with
```
newparams
```

The full configuration for a single geodesic is stored in:

```@docs
CarterGeodesicParams
```