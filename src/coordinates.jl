
"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Sigma = r^2 + a^2 \\cos^2( \\theta ).
```
"""
@vec_eq Σ(r, a, θ) = r^2 + (a * cos(θ))^2


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.3):

```math
\\Delta = r^2 - 2 M r + a^2.
```
"""
@vec_eq Δ(M, r, a) = r^2 - 2 * M * r + a^2


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.3):

```math
A = (r^2 + a^2)^2 - a^2 \\Delta \\sin^2 ( \\theta ).
```
"""
@vec_eq A(M, r, a, θ) = (r^2 + a^2)^2 - a^2 * Δ(M, r, a) * sin(θ)^2


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.10):

```math
T = E \\left( r^2 + a^2 \\right) - L * a.
```
"""
@vec_eq T(E, L, r, a) = E * (r^2 + a^2) - L * a


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_r = T^2 - \\Delta \\left[ (L - a E)^2 + Q \\right]
```
"""
@vec_eq Vr(E, L, M, Q, r, a) = T(E, L, r, a)^2 - Δ(M, r, a) * ((L - a * E)^2 + Q)


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_\\theta = 
    Q + \\cos^2 (\\theta) \\left[ a^2 E^2 - \\frac{L^2}{\\sin^2 (\\theta) } \\right].
```
"""
@vec_eq Vθ(E, L, Q, a, θ) = Q + cos(θ)^2 * ((a * E)^2 - (L / sin(θ))^2)


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_1 = 
    1 + \\sqrt[3]{1 - \\frac{a^2}{M^2}} 
    \\left[ 
        \\sqrt[3]{1 + \\frac{a}{M}} + \\sqrt[3]{1 - \\frac{a}{M}} 
    \\right].
```
"""
@vec_eq Z₁(M, a) = 1 + ∛(1 - (a / M)^2) * (∛(1 + (a / M)) + ∛(1 - (a / M)))


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.21):

```math
Z_2 = \\sqrt{\\frac{3a^2}{M^2} + Z_1^2}.
```
"""
@vec_eq Z₂(M, a) = √(3(a / M)^2 + Z₁(M, a)^2)


"""
    $(TYPEDSIGNATURES)

Radius of marginally stable orbit.

From Bardeen et al. (1972) eq. (2.21):

```math
r_\\text{ms} = M \\left\\{ 3 + Z_2 \\pm \\sqrt{(3 - Z_1)(3 + Z_1 + 2 Z_2)} \\right\\}.
```

Depending on whether the [`AbstractSpinDirection`](@ref) is `ContraRotating` or 
`CoRotating`, changes the sign in the equation above.
"""
@inline rms(M, a; ± = +) =
    M * (3 + Z₂(M, a) ± √((3 - Z₁(M, a)) * (3 + Z₁(M, a) + 2 * Z₂(M, a))))
@inline rms(p::GeodesicParams) = p.rms
@inline rms(s::BHSetup) = s.a < 0.0 ? rms(s.M, s.a) : rms(s.M, s.a; ± = -)


"""
    $(TYPEDSIGNATURES)

Event horizon radius for a black hole of mass ``M`` and spin ``a``.

From Bardeen et al. (1972) eq. (2.7):

```math
R_0 = M + \\sqrt{M^2 + a^2}.
```
"""
R₀(M, a) = M + √(M^2 - a^2)
R₀(s::BHSetup) = R₀(s.M, s.a)

"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

The ``t`` compontent of the equation of motion for a photon around a black hole, multiplied 
by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9d):

```math
\\Sigma \\frac{\\text{d}t}{\\text{d}\\lambda} =
    - a \\left( a E \\sin^2 \\theta - L \\right)
    + \\frac{\\left( r^2 + a^2 \\right) T}{\\Delta}.
```
"""
@vec_eq function Σδt_δλ(E, L, M, r, a, θ)
    -a * (a * E * sin(θ)^2 - L) + (r^2 + a^2) * T(E, L, r, a) / Δ(M, r, a)
end


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

The ``r`` compontent of the equation of motion for a photon around a black hole, multiplied 
by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9a):

```math
\\Sigma \\frac{\\text{d}r}{\\text{d}\\lambda} =
    \\pm \\sqrt{\\lvert V_r \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied 
in [`δ`](@ref).
"""
@vec_eq function Σδr_δλ(E, L, M, Q, r, a)
    V = Vr(E, L, M, Q, r, a)
    √(V * sign(V))
end


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

The ``\\theta`` compontent of the equation of motion for a photon around a black hole, 
multiplied by ``\\Sigma``.

Modified from Bardeen et al. (1972) eq. (2.9b):

```math
\\Sigma \\frac{\\text{d}\\theta}{\\text{d}\\lambda} =
    \\pm \\sqrt{\\lvert V_\\theta \\rvert},
```

where, for implementation reason, the sign is always positive. Instead, the sign is applied 
in [`δ`](@ref).
"""
@vec_eq function Σδθ_δλ(E, L, Q, a, θ)
    V = Vθ(E, L, Q, a, θ)
    √(V * sign(V))
end


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

The ``\\phi`` compontent of the equation of motion for a photon around a black hole, 
multiplied by ``\\Sigma``.

From Bardeen et al. (1972) eq. (2.9c):

```math
\\Sigma \\frac{\\text{d}\\phi}{\\text{d}\\lambda} =
    - \\frac{L}{\\sin^2 \\theta} - aE
    + \\frac{aT}{\\Delta}.
```
"""
@vec_eq function Σδϕ_δλ(E, L, M, r, a, θ)
    (L / sin(θ)^2) - (a * E) + a * T(E, L, r, a) / Δ(M, r, a)
end


"""
    $(TYPEDSIGNATURES)

Derivative of `x`, a [`FourVector`](@ref), with respect to ``\\lambda``, for the Kerr 
spacetime around a black hole described by `p`, of type [`GeodesicParams`](@ref).

The parameter `signs` is a mutable structure with two components, where the first is the 
sign of ``V_r``, and the second ``V_\\theta``.

Returns a [`FourVector`](@ref) with components

```math
\\left(
    \\frac{\\text{d}t}{\\text{d}\\lambda},
    \\frac{\\text{d}r}{\\text{d}\\lambda},
    \\frac{\\text{d}\\theta}{\\text{d}\\lambda}, 
    \\frac{\\text{d}\\phi}{\\text{d}\\lambda}
\\right).
```
"""
@inline function δ(x::T, p::GeodesicParams)::T where {T<:AbstractVector}
    Σ₀ = Σ(x, p)
    T(
        Σδt_δλ(x, p) / Σ₀,
        p.r_sign * Σδr_δλ(x, p) / Σ₀,
        p.θ_sign * Σδθ_δλ(x, p) / Σ₀,
        Σδϕ_δλ(x, p) / Σ₀
    )
end

"""
    $(TYPEDSIGNATURES)

Inplace variant of [`δ`](@ref).
"""
@inline function δ!(res, x, p::GeodesicParams)
    Σ₀ = Σ(x, p)
    @inbounds begin
        res[1] = Σδt_δλ(x, p) / Σ₀
        res[2] = p.r_sign * Σδr_δλ(x, p) / Σ₀
        res[3] = p.θ_sign * Σδθ_δλ(x, p) / Σ₀
        res[4] = Σδϕ_δλ(x, p) / Σ₀
    end
end
