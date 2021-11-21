
"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Sigma = r^2 + a^2 \\cos^2( \\theta ).
```
"""
Σ(r, a, θ) = r^2 + (a * cos(θ))^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Delta = r^2 - 2 M r + a^2.
```
"""
Δ(M, r, a) = r^2 - 2 * M * r + a^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.3):

```math
A = (r^2 + a^2)^2 - a^2 \\Delta \\sin^2 ( \\theta ).
```
"""
A(M, r, a, θ) = (r^2 + a^2)^2 - a^2 * Δ(M, r, a) * sin(θ)^2


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_1 = 
    1 + \\sqrt[3]{1 - \\frac{a^2}{M^2}} 
    \\left[ 
        \\sqrt[3]{1 + \\frac{a}{M}} + \\sqrt[3]{1 - \\frac{a}{M}} 
    \\right].
```
"""
Z₁(M, a) = 1 + ∛(1 - (a / M)^2) * (∛(1 + (a / M)) + ∛(1 - (a / M)))


"""
    $(TYPEDSIGNATURES)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_2 = \\sqrt{\\frac{3a^2}{M^2} + Z_1^2}.
```
"""
Z₂(M, a) = √(3(a / M)^2 + Z₁(M, a)^2)


"""
    $(TYPEDSIGNATURES)

Radius of marginally stable orbit.

From Bardeen et al. (1972) eq. (2.21):

```math
r_\\text{ms} = M \\left\\{ 3 + Z_2 \\pm \\sqrt{(3 - Z_1)(3 + Z_1 + 2 Z_2)} \\right\\}.
```

The choice of ``\\pm`` is chosen by the sign of ``a``.
"""
@inline function rms(M, a, ±)
    M * (3 + Z₂(M, a) ± √((3 - Z₁(M, a)) * (3 + Z₁(M, a) + 2 * Z₂(M, a))))
end
@inline function rms(M, a) 
    a > 0.0 ? rms(M, a, -) : rms(M, a, +)
end 
@inline function rms(s::BHSetup)
    rms(s.metric.M, s.metric.a)
end

"""
    $(TYPEDSIGNATURES)

Event horizon radius for a black hole of mass ``M`` and spin ``a``.

From Bardeen et al. (1972) eq. (2.7):

```math
R_0 = M + \\sqrt{M^2 + a^2}.
```
"""
@inline R₀(M, a) = M + √(M^2 - a^2)
@inline R₀(s::BHSetup) = R₀(s.metric.M, s.metric.a)
