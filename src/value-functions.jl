
struct ValueFunction{F}
    func::F
    initval::Float64

    ValueFunction(func::F, initval) where {F} = new{F}(func, initval)
    ValueFunction(func) = ValueFunction(func, NaN64)
end

@inline function (vf::ValueFunction)(last_λ, last_u, last_v, p, s::BHSetup, d::AccretionDisk)::Float64 where {T}
    @inbounds if last_λ < s.λhigh && d.r_inner < last_u[2] < d.r_outer
        return vf.func(0.0, last_λ, last_u, last_v, p, d)
    end
    vf.initval
end

@inline function (vf::ValueFunction)(sol, s::BHSetup{CarterBoyerLindquist{T}}, d::AccretionDisk)::Float64 where {T}
    @inbounds vf(sol.t[end], sol.u[end], nothing, sol.prob.p, s, d)
end
@inline function (vf::ValueFunction)(sol, s, d::AccretionDisk)::Float64
    @inbounds vf(sol.t[end], sol.u[end].x[2], sol.u[end].x[1], sol.prob.p, s, d)
end

# no disk specialisations
@inline function (vf::ValueFunction)(sol, s::BHSetup{CarterBoyerLindquist{T}}, d::Nothing)::Float64 where {T}
    vf.func(vf.initval, sol.t[end], sol.u[end], nothing, sol.prob.p, d)
end
@inline function (vf::ValueFunction)(sol, s, d::Nothing)::Float64
    vf.func(vf.initval, sol.t[end], sol.u[end].x[2], sol.u[end].x[1], sol.prob.p, d)
end

@inline function Base.:∘(vf1::ValueFunction, vf2::ValueFunction)
    ValueFunction(
        (val, last_λ, last_u, last_v, p, d) ->
            vf1.func(vf2.func(val, last_λ, last_u, last_v, p, d), last_λ, last_u, last_v, p, d)
    )
end

# a few commonly used value functions

const geometry = ValueFunction((val, λ, u, v, p, d) -> 1.0)

export ValueFunction, geometry, redshift
