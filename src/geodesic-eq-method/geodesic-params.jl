
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
@with_kw struct GeodesicParams{M,V,T} <: AbstractGeodesicParams{V,T}
    @deftype T
    metric::M

    "Initial radial velocity angle."
    ϕv₀
    "Initial azimuthal velocity angle."
    θv₀

    chart_inner_radius::Float32

    "Storage parameter, used when custom storage for callbacks is required."
    storage::V = nothing
end

"""
    $(TYPEDSIGNATURES)
"""
GeodesicParams(α, β, s::BHSetup) = GeodesicParams(α, β, s::BHSetup, nothing)
function GeodesicParams(α, β, s::BHSetup{M}, storage) where {M}
    GeodesicParams{M}(
        metric = s.metric,
        # TODO: impact parameter scaling
        # applied here is just a heuristic, noticing that scaling with s.r₀^2 approximately reproduces
        # the impact paramter mapping of Bardeen and Cunningham
        ϕv₀ = α / s.r₀^2,
        θv₀ = β / s.r₀^2,
        chart_inner_radius = R₀(metric.M, metric.a),
        storage = storage
    )
end

"""
    $(TYPEDSIGNATURES)
"""
function newparams(p::GeodesicParams, θ, r, α, β, δα)::GeodesicParams
    @set(p.ϕv₀ = α + δα / s.r₀^2)
end
