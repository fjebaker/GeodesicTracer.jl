"""
    $(TYPEDEF)

The struct has the following fields all of type `T`:
$(FIELDS)

This struct is not mean to be instatiated by the user, and is used by the integrator as the 
integration parameter; it is immutable, and if value changes are desired, may be altered 
with `@set` from [Accessors.jl](https://juliaobjects.github.io/Accessors.jl/stable/).

Note it is currently possible for ``\\lvert a \\rvert > 1``, and for the `rms` to be
set incorreclty for a `a`. Doing this is not advised, and generally, instances of this
type should not be created or modified by the user.
"""

@with_kw struct GeodesicParams{V,T}
    @deftype T
    M
    E
    a

    L
    Q
    rms::Float32
    R₀::Float32

    "Used to track the effective potential ``V_r`` sign."
    r_sign::Float32 = -1.0
    "Used to track the effective potential ``V_\\theta`` sign."
    θ_sign::Float32 = 0.0

    "Proper time parameter when ``V_r`` sign changes (set by the integrator callback)."
    λr_change = 0.0
    "Proper time parameter when ``V_\\theta`` sign changes (set by the integrator callback)."
    λθ_change = 0.0

    "Storage parameter, used when custom storage for callbacks is required."
    storage::V = nothing
end

# forward declare 
function LQ end

"""
    $(TYPEDSIGNATURES)
"""
GeodesicParams(α, β, s::BHSetup) = GeodesicParams(α, β, s::BHSetup, nothing)
function GeodesicParams(α, β, s::BHSetup, storage)
    l, q = LQ(s.M, s.r₀, s.a, s.θ₀, α, β)
    GeodesicParams(
        M = s.M,
        E = s.E,
        a = s.a,
        L = l,
        Q = q,
        rms = rms(s),
        R₀ = R₀(s.M, s.a),
        θ_sign = β > 0 ? -1.0 : 1.0,
        storage = storage
    )
end

"""
    $(TYPEDSIGNATURES)

Return a new [`GeodesicParams`](@ref) with the `r_sign` parameter sign flipped, along with the proper time `λ` stored in `λr_change`.
"""
flip_rsign(λ, p::GeodesicParams) = @set(@set(p.r_sign = -p.r_sign).λr_change = λ)

"""
    $(TYPEDSIGNATURES)

Return a new [`GeodesicParams`](@ref) with the `θ_sign` parameter sign flipped, along with the proper time `λ` stored in `λθ_change`.
"""
flip_θsign(λ, p::GeodesicParams) = @set(@set(p.θ_sign = -p.θ_sign).λθ_change = λ)

"""
    $(TYPEDSIGNATURES)

Casts [`GeodesicParams`](@ref) of type `{S,T}` to `{S,NewType}`.
"""
changetype(NewType::Type, p::GeodesicParams) = GeodesicParams(;
    (f => NewType(getfield(p, f)) for f ∈ fieldnames(GeodesicParams) if f != :storage)...
)
1.0
