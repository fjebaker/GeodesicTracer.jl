"""
    $(TYPEDEF)

The struct has the following fields all of type `T`:
$(FIELDS)

A Carter-Boyer-Lindquist method specialisation of [`GeodesicParams`](@ref).
"""

@with_kw struct CarterGeodesicParams{V,T} <: AbstractGeodesicParams{V,T}
    @deftype T
    metric::CarterBoyerLindquist{T}

    L
    Q
    rms::Float32

    chart_inner_radius::Float32

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

CarterGeodesicParams(α, β, s::BHSetup) = CarterGeodesicParams(α, β, s::BHSetup, nothing)
function CarterGeodesicParams(α, β, s::BHSetup, storage)
    metric = s.metric
    l, q = LQ(metric.M, s.r₀, metric.a, s.θ₀, α, β)
    CarterGeodesicParams(
        metric = metric,
        L = l,
        Q = q,
        rms = rms(s),
        chart_inner_radius = R₀(metric.M, metric.a),
        θ_sign = β > 0 ? -1.0 : 1.0,
        storage = storage
    )
end

"""
    $(TYPEDSIGNATURES)

Return a new [`CarterGeodesicParams`](@ref) with the `r_sign` parameter sign flipped, along with the proper time `λ` stored in `λr_change`.
"""
flip_rsign(λ, p::CarterGeodesicParams) = @set(@set(p.r_sign = -p.r_sign).λr_change = λ)

"""
    $(TYPEDSIGNATURES)

Return a new [`CarterGeodesicParams`](@ref) with the `θ_sign` parameter sign flipped, along with the proper time `λ` stored in `λθ_change`.
"""
flip_θsign(λ, p::CarterGeodesicParams) = @set(@set(p.θ_sign = -p.θ_sign).λθ_change = λ)


"""
    $(TYPEDSIGNATURES)

Carter method specialisation.
"""
function makeprobfunc(s::BHSetup{CarterBoyerLindquist{T}}, α_range, β, num) where {T}
    α = α_range[1]
    δα = (α_range[2] - α) / num
    metric = s.metric

    (prob, i, repeat) -> begin
        p = prob.p
        r = prob.u0[2]
        θ = prob.u0[3]

        l, q = LQ(metric.M, r, metric.a, θ, α + i * δα, β)
        p = @set(@set(p.L = l).Q = q)
        remake(prob, p = p)
    end
end
