# Coordinate functions

```@meta
CurrentModule = GeodesicTracer
```

## Coordinate functions

The Kerr metric in the Boyer-Lindquist coordinates is written

```math
\begin{aligned}
\text{d}s^2 = & - \left( 1 - \frac{2Mr}{\Sigma} \right) \text{d}t^2 \\
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
A
```

The event horizon and marginally stable / innermost stable circular orbit (ISCO) is given by: 

```@docs
R₀
rms
Z₁
Z₂
```

### Equation of motion for null geodesics

For a photon travelling along null geodesics, the governing components are:

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

### Constants of motion

The constants of motion are:

```@docs
T
S
LQ
```

### Convenience functions 

```@docs
δ
δ!
```

Where the observer angles for impact parameters $\alpha$ and $\beta$ 

```@docs
sinΦsinΨ
```