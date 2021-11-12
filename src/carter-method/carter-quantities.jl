"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.10):

```math
T = E \\left( r^2 + a^2 \\right) - L * a.
```
"""
T(E, L, r, a) = E * (r^2 + a^2) - L * a


"""
    $(TYPEDSIGNATURES)
    $(FUNCTIONNAME)(x::FourVector, p::GeodesicParams)
    $(FUNCTIONNAME)(x::AbstractVector, p::GeodesicParams)

From Bardeen et al. (1972) eq. (2.10), for the special case of a null-geodesic ``\\mu = 0``:

```math
V_r = T^2 - \\Delta \\left[ (L - a E)^2 + Q \\right]
```
"""
Vr(E, L, M, Q, r, a) = T(E, L, r, a)^2 - Δ(M, r, a) * ((L - a * E)^2 + Q)


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
Vθ(E, L, Q, a, θ) = Q + cos(θ)^2 * ((a * E)^2 - (L / sin(θ))^2)


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
function Σδt_δλ(E, L, M, r, a, θ)
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
function Σδr_δλ(E, L, M, Q, r, a)
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
function Σδθ_δλ(E, L, Q, a, θ)
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
function Σδϕ_δλ(E, L, M, r, a, θ)
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
@inline function δ(x::T, p)::T where {T<:AbstractVector}
    metric = p.metric
    Σ₀ = Σ(x[2], metric.a, x[3])
    T(
        Σδt_δλ(metric.E, p.L, metric.M, x[2], metric.a, x[3]) / Σ₀,
        p.r_sign * Σδr_δλ(metric.E, p.L, metric.M, p.Q, x[2], metric.a) / Σ₀,
        p.θ_sign * Σδθ_δλ(metric.E, p.L, p.Q, metric.a, x[3]) / Σ₀,
        Σδϕ_δλ(metric.E, p.L, metric.M, x[2], metric.a, x[3]) / Σ₀
    )
end

"""
    $(TYPEDSIGNATURES)

Inplace variant of [`δ`](@ref).
"""
@inline function δ!(res, x)
    metric = p.metric
    Σ₀ = Σ(x[2], metric.a, x[3])
    @inbounds begin
        res[1] = Σδt_δλ(metric.E, p.L, metric.M, x[2], metric.a, x[3]) / Σ₀
        res[2] = p.r_sign * Σδr_δλ(metric.E, p.L, metric.M, p.Q, x[2], metric.a) / Σ₀
        res[3] = p.θ_sign * Σδθ_δλ(metric.E, p.L, p.Q, metric.a, x[3]) / Σ₀
        res[4] = Σδϕ_δλ(metric.E, p.L, metric.M, x[2], metric.a, x[3]) / Σ₀
    end
end


"""
Carter method specialsation.
"""
function rayintegrator(u, p::CarterGeodesicParams, λ)
    δ(u, p)
end

"""
Carter method specialsation.
"""
function rayintegrator!(du, u, p::CarterGeodesicParams, λ)
    δ!(du, u, p)
end
function rayintegrator!(du, u, p::AbstractArray{CarterGeodesicParams}, λ)
    # gpu version requires parameter unpacking
    δ!(du, u, p[1])
end
