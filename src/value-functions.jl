
struct ValueFunction{F}
    func::F
end

@inline function (vf::ValueFunction)(sol, s::BHSetup, d::AccretionDisk)::Float64
    last_λ = sol.t[end]
    last_u = sol.u[end]
    p = sol.prob.p
    @inbounds if last_λ < s.λhigh && d.r_inner < last_u[2] < d.r_outer
        return vf.func(0.0, last_λ, last_u, p, d)
    end
    0.0
end

# no disk specialisation
@inline function (vf::ValueFunction)(sol, s::BHSetup, d::Nothing)::Float64
    last_λ = sol.t[end]
    last_u = sol.u[end]
    p = sol.prob.p
    vf.func(0.0, last_λ, last_u, p, d)
end


@inline function Base.:∘(vf1::ValueFunction, vf2::ValueFunction)
    ValueFunction(
        (val, last_λ, last_u, p::GeodesicParams, d) ->
            vf1.func(vf2.func(val, last_λ, last_u, p, d), last_λ, last_u, p, d)
    )
end

# a few commonly used value functions

const redshift = ValueFunction(
    (val, λ, u, p::GeodesicParams, d) -> begin
        @inbounds if u[2] > p.rms
            return reg_pdotu_inv(u, p)
        else
            return plg_pdotu_inv(u, p, last_t < p.λr_change ? -1 : 1)
        end
    end
)

const geometry = ValueFunction((val, λ, u, p::GeodesicParams, d) -> 1.0)


export ValueFunction, geometry, redshift
