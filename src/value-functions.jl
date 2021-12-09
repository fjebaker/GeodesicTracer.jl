
struct ValueFunction{F}
    func::F
end

@inline function (vf::ValueFunction)(last_λ, last_u, last_v, p, s::BHSetup, d::AccretionDisk)::Float64 where {T}
    @inbounds if last_λ < s.λhigh && d.r_inner < last_u[2] < d.r_outer
        return vf.func(0.0, last_λ, last_u, last_v, p, d)
    end
    0.0
end

@inline function (vf::ValueFunction)(sol, s::BHSetup{CarterBoyerLindquist{T}}, d::AccretionDisk)::Float64 where {T}
    @inbounds vf(sol.t[end], sol.u[end], nothing, sol.prob.p, s, d)
end
@inline function (vf::ValueFunction)(sol, s, d::AccretionDisk)::Float64
    @inbounds vf(sol.t[end], sol.u[end].x[2], sol.u[end].x[1], sol.prob.p, s, d)
end

# no disk specialisations
@inline function (vf::ValueFunction)(sol, s::BHSetup{CarterBoyerLindquist{T}}, d::Nothing)::Float64 where {T}
    vf.func(0.0, sol.t[end], sol.u[end], nothing, sol.prob.p, d)
end
@inline function (vf::ValueFunction)(sol, s, d::Nothing)::Float64
    vf.func(0.0, sol.t[end], sol.u[end].x[2], sol.u[end].x[1], sol.prob.p, d)
end

@inline function Base.:∘(vf1::ValueFunction, vf2::ValueFunction)
    ValueFunction(
        (val, last_λ, last_u, last_v, p, d) ->
            vf1.func(vf2.func(val, last_λ, last_u, last_v, p, d), last_λ, last_u, last_v, p, d)
    )
end

# a few commonly used value functions

@inline function redshift_function(val, λ, u, v, p::CarterGeodesicParams, d)
    # this function needs a specialsation for GeodesicParams
    # since at the moment reg_pdotu_inv and plg_pdotu_inv are Carter specific
    @inbounds if u[2] > rms(p.metric.M, p.metric.a)
        return reg_pdotu_inv(u, p)
    else
        return plg_pdotu_inv(u, p, λ < p.λr_change ? -1 : 1)
    end
end

@inline function redshift_function(val, λ, u, v, p::GeodesicParams, d)
    # this function needs a specialsation for GeodesicParams
    # since at the moment reg_pdotu_inv and plg_pdotu_inv are Carter specific
    @inbounds if u[6] > rms(p.metric.M, p.metric.a)
        return reg_pdotu_inv(u, p)
    else
        return plg_pdotu_inv(u, p, λ < p.λr_change ? -1 : 1)
    end
end

const redshift = ValueFunction(
    redshift_function
)

const geometry = ValueFunction((val, λ, u, v, p, d) -> 1.0)

export ValueFunction, geometry, redshift
