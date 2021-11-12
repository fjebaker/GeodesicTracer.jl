
"""
    $(TYPEDSIGNATURES)

Modified from Cunningham et al. (1975) eq. (A2a):

```math
e^\\nu = \\sqrt{\\frac{\\Delta \\Sigma}{A}}.
```
"""
eⱽ(M, r, a, θ) = √(Σ(r, a, θ) * Δ(M, r, a) / A(M, r, a, θ))


"""
    $(TYPEDSIGNATURES)

Modified from Cunningham et al. (1975) eq. (A2b):

```math
e^\\Phi = \\sin \\theta \\sqrt{\\frac{A}{\\Sigma}}.
```
"""
eᶲ(M, r, a, θ) = sin(θ) * √(A(M, r, a, θ) / Σ(r, a, θ))


"""
    $(TYPEDSIGNATURES)

From Cunningham et al. (1975) eq. (A2c):

```math
\\omega = \\frac{2 a M r}{A}.
```
"""
ω(M, r, a, θ) = 2 * a * M * r / A(M, r, a, θ)


"""
    $(TYPEDSIGNATURES)

Coordinate angular velocity of an accreting gas. 

Taken from Cunningham et al. (1975) eq. (A7b):

```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} + r_e^{3/2}}.
```

# Notes

Fanton et al. (1997) use
```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} \\pm r_e^{3/2}},
```
where the sign is dependent on co- or contra-rotation. This function may be extended in the future to support this definition.
"""
Ωₑ(M, r, a) = √M / (r^1.5 + a * √M)


"""
    $(TYPEDSIGNATURES)

Velocity of an accreting gas in a locally non-rotating reference frame (see Bardeen et al. 1973).

Taken from Cunningham et al. (1975) eq. (A7b):

```math
V_e = (\\Omega_e - \\omega) e^{\\Phi - \\nu}.
```
"""
Vₑ(M, r, a, θ) = (Ωₑ(M, r, a) - ω(M, r, a, θ)) * eᶲ(M, r, a, θ) / eⱽ(M, r, a, θ)


"""
    $(TYPEDSIGNATURES)

Angular momentum of an accreting gas within ``r_ms``.

Taken from Cunningham et al. (1975) eq. (A11b):

```math
L_e = \\sqrt{M} \\frac{
        r_{\\text{ms}}^2 - 2 a \\sqrt{M r_{\\text{ms}}} + a^2
    }{
        r_{\\text{ms}}^{3/2} - 2 M \\sqrt{r_{\\text{ms}}} + a \\sqrt{M}
    }.
```
"""
Lₑ(M, rms, a) = √M * (rms^2 - 2 * a * √(M * rms) + a^2) / (rms^1.5 - 2 * M * √rms + a * √M)


"""
    $(TYPEDSIGNATURES)

Taken from Cunningham et al. (1975) eq. (A12e):

```math
H = \\frac{2 M r_e - a \\lambda_e}{\\Delta},
```
where we distinguing ``r_e`` as the position of the accreting gas. 
"""
H(M, rms, r, a) = (2 * M * r - a * Lₑ(M, rms, a)) / Δ(M, r, a)


"""
    $(TYPEDSIGNATURES)

Taken from Cunningham et al. (1975) eq. (A11c):

```math
\\gamma_e = \\sqrt{1 - \\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }}.
```
"""
γₑ(M, rms) = √(1 - (2 * M) / (3 * rms))


"""
    $(TYPEDSIGNATURES)

Taken from Cunningham et al. (1975) eq. (A12b):

```math
u^r = - \\sqrt{\\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }} \\left(
        \\frac{ r_{\\text{ms}} }{r_e} - 1
    \\right)^{3/2}.
```
"""
uʳ(M, rms, r) = -√((2 * M) / (3 * rms)) * (rms / r - 1)^1.5


"""
    $(TYPEDSIGNATURES)

Taken from Cunningham et al. (1975) eq. (A12c):

```math
u^\\phi = \\frac{\\gamma_e}{r_e^2} \\left( 
        L_e + aH 
    \\right).
```
"""
uᶲ(M, rms, r, a) = γₑ(M, rms) / r^2 * (Lₑ(M, rms, a) + a * H(M, rms, r, a))


"""
    $(TYPEDSIGNATURES)

Taken from Cunningham et al. (1975) eq. (A12b):

```math
u^t = \\gamma_e \\left(
        1 + \\frac{2 M (1 + H)}{r_e}
    \\right).
```
"""
uᵗ(M, rms, r, a) = γₑ(M, rms) * (1 + 2 * M * (1 + H(M, rms, r, a)) / r)


"""
    $(TYPEDSIGNATURES)

BLANK
"""
function reg_pdotu_inv(L, M, r, a, θ)
    (eⱽ(M, r, a, θ) * √(1 - Vₑ(M, r, a, θ)^2)) / (1 - L * Ωₑ(M, r, a))
end

"""
    $(TYPEDSIGNATURES)

BLANK
"""
function plg_pdotu_inv(E, L, M, Q, rms, r, a, sign_r)
    inv(
        uᵗ(M, rms, r, a) - uᶲ(M, rms, r, a) * L -
        sign_r * uʳ(M, rms, r) * Σδr_δλ(E, L, M, Q, r, a) / Δ(M, r, a)
    )
end