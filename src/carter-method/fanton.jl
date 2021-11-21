
"""
    $(TYPEDSIGNATURES)

From Fanton et al. (1997), eq. (76):

```math
S = 1 + \\alpha \\omega \\sin \\theta.
```
"""
S(θ, α, ω) = 1 + α * ω * sin(θ)


"""
    $(TYPEDSIGNATURES)

Calculates and returns the observer's angles ``\\sin \\Theta`` and ``\\sin \\Phi``, where 
the parameters in the function signature correspond to the Bardeen et al. (1972) paper.

From Fanton et al. (1997), eq. (74) and eq. (75):

```math
\\begin{align*}
\\sin \\Theta &=
    \\frac{\\alpha \\Sigma}{\\sqrt{A}} 
    \\left\\{ 
            \\beta^2 + \\left( \\alpha + a \\sin \\theta \\right)^2
            + \\frac{
                A S^2 - \\left( r^2 + a^2 + a \\alpha \\sin \\theta \\right) 
            }{\\Delta} 
    \\right\\}, \\\\
\\sin \\Phi &= 
    -\\frac{\\alpha \\Sigma \\sqrt{\\Delta}}{A S \\sin\\Theta}.
\\end{align*}
```
"""
function sinΦsinΨ(Σ₀, sinθ, A₀, Δ₀, S₀, r, a, α, β)
    # calc 1
    sinΦ =
        ((α * Σ₀) / √A₀) /
        sqrt(β^2 + (α + a * sinθ)^2 + (A₀ * S₀^2 - (r^2 + a^2 + a * α * sinθ)^2) / Δ₀)

    # calc 2
    sinΨ = -(α * Σ₀ * √Δ₀) / (S₀ * A₀ * sinΦ)

    # return
    (sinΦ, sinΨ)
end


"""
    $(TYPEDSIGNATURES)

Calculates conserved quantities
    - angular momentum ``L`` 
    - Carter parameter ``Q``

for a photon described with position described by `x` in a Kerr spacetime given by 
`p`.

From Fanton et al. (1997), eq. (69):

```math
L = \\frac{\\Upsilon_1}{\\Upsilon_2},
```

where

```math
\\begin{align*}
    \\Upsilon_1 &= \\sin \\theta \\sin \\Phi \\sin \\Theta, \\\\
    \\Upsilon_2 &= \\frac{\\Sigma \\sqrt{\\Delta}}{A} + \\omega \\Upsilon_1,
    \\omega = \\frac{2 a r}{A},
\\end{align*}
```

taken from eq. (72) and eq. (73).

From Fanton et al. (1997), eq. (70):

```math
Q = 
    \\frac{P^2}{\\Delta} 
    - \\left( \\lambda - a \\right)^2 
    - \\frac{\\Sigma^2}{A} \\left( \\frac{\\cos \\Phi}{\\Upsilon_2} \\right)^2,
```

and

```math
P = \\left( r^2 + a^2 - a L \\right),
```

taken from eq. (71).
"""
function LQ(M, r, a, θ, α, β)
    # value cache
    Σ₀ = Σ(r, a, θ)
    sinθ = sin(θ)
    A₀ = A(M, r, a, θ)
    Δ₀ = Δ(M, r, a)
    ω = 2.0 * a * r / A₀
    S₀ = S(θ, α, ω)

    # calculations
    sinΦ, sinΨ = sinΦsinΨ(Σ₀, sinθ, A₀, Δ₀, S₀, r, a, α, β)
    Υ₁ = sinθ * sinΦ * sinΨ
    Υ₂ = (Σ₀ * √Δ₀ / A₀) + ω * Υ₁

    L = Υ₁ / Υ₂
    P = (r^2 + a^2 - (a * L))
    Q = (P^2 / Δ₀) - (L - a)^2 - ((Σ₀^2 / A₀) * (cos(asin(sinΨ)) / Υ₂)^2)

    (L, Q)
end
