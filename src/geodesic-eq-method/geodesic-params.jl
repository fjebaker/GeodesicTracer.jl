using LinearAlgebra
using RecursiveArrayTools

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

function gθθ(r, θ, M, a)
    a^2 * cos(θ)^2 + r^2
end

function gtt(r, θ, M, a)
    -(2 * 2*M * r * a * sin(θ)^2) / r^2 + a^2 * cos(θ)^2
end

function gϕϕ(r, θ, M, a)
    sin(θ)^2 * (r^2 + a^2 + sin(θ)^2 * 2*M * r * a^2/(r^2 + a^2 * cos(θ)^2))
end


    

"""
    $(TYPEDSIGNATURES)
"""
GeodesicParams(α, β, s::BHSetup) = GeodesicParams(α, β, s::BHSetup, nothing)
function GeodesicParams(α, β, s::BHSetup{M}, storage) where {M}
    GeodesicParams(
        metric = s.metric,
        # TODO: impact parameter scaling
        # applied here is just a heuristic, noticing that scaling with s.r₀^2 approximately reproduces
        # the impact paramter mapping of Bardeen and Cunningham
        ϕv₀ = -α / gϕϕ(s.r₀, s.θ₀, s.metric.M, s.metric.a),
        θv₀ = -β / gθθ(s.r₀, s.θ₀, s.metric.M, s.metric.a),
        # scale inner chart, since otherwise infinite affine parameter as approaching event horizon
        chart_inner_radius = 1.06 * R₀(s.metric.M, s.metric.a),
        storage = storage
    )
end


"""
    $(TYPEDSIGNATURES)
"""
function makeprobfunc(s::BHSetup, α_range, β, num)
    α = α_range[1]
    δα = (α_range[2] - α) / num
    metric = s.metric
    
    (prob, i, repeat) -> begin
        p = prob.p
        x = prob.u0.x[1]
        
        p = @set(p.ϕv₀ = -(α + i*δα) / gϕϕ(x[2], x[3], metric.M, metric.a))

        # calculate new velocity vector
        vtemp = (0.0, -1.0, p.θv₀, p.ϕv₀)
        v = SVector(null_constrain(x, vtemp, metric), -1.0, p.θv₀, p.ϕv₀)

        remake(prob, p = p, u0 = ArrayPartition(v, x))
    end
end