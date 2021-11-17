"""
    $(TYPEDSIGNATURES)

Integration problem. Here, `du` represents velocity 4-vector

```math
\\left(
    \\frac{\\text{d}t}{\\text{d}\\lambda},
    \\frac{\\text{d}r}{\\text{d}\\lambda},
    \\frac{\\text{d}\\theta}{\\text{d}\\lambda},
    \\frac{\\text{d}\\phi}{\\text{d}\\lambda}
\\right),
```

and `u` is the regular 4-vector

```math
\\left( t, r, \\theta, \\phi \\right).
```

These are calculated by [`δ`](@ref).

# Developer notes

This function can be dispatched over `p` if needed.
"""
function rayintegrator(v, u, p::GeodesicParams, λ)
    SVector(geodesic_eq(u, v, p.metric)...)
end

"""
    $(TYPEDSIGNATURES)

Inplace variant of [`rayintegrator`](@ref).
"""
function rayintegrator!(du, u, p::GeodesicParams, λ)
    x = @inbounds @view(u[1:4])
    v = @inbounds @view(u[5:8])
    du[1:4] .= v
    du[5:8] .= geodesic_eq(x, v, p.metric)
end
rayintegrator!(du, u, p::AbstractArray{GeodesicParams}, λ) = rayintegrator!(du, u, p[1], λ)


"""
    $(TYPEDSIGNATURES)

Integrate a geodesic with impact parameters ``\\alpha``, ``\\beta``, for a Kerr spacetime 
described by `init_p`, with geometric setup `s`. 

Saves the entire geodesic path, unless `save_geodesics=false`. If `disk` is supplied,
will create a callback variant with [`wrapcallback`](@ref), which checks whether the
geodesic intersects the disk.
"""
function calcgeodesic(
    α::AbstractFloat,
    β::AbstractFloat,
    s::BHSetup;
    save_geodesics = true,
    disk = nothing
)
    # single geodesic method
    integrategeodesic(
        s,
        SingleParams(
            α = α,
            β = β,
            save_geodesics = save_geodesics,
            callback = wrapcallback(s, disk)
        ),
        storage = isnothing(disk) ? nothing : 0.0
    )
end


"""
    $(TYPEDSIGNATURES)

Specialisation for ensemble calculations across a range of `num` impact parameters 
``\\alpha`` between `α_range[1]` and `α_range[2]`.
"""
function calcgeodesic(
    α_range::Tuple{Float64,Float64},
    num,
    β::AbstractFloat,
    s::BHSetup;
    save_geodesics = true,
    disk = nothing
)
    
    integrategeodesic(
        s,
        ParallelParams(
            trajectories = num,
            β = β,
            α = α_range[1],
            α_end = α_range[2],
            save_geodesics = save_geodesics,
            ensemble = EnsembleThreads(),
            probfunc = makeprobfunc(s, α_range, β, num),
            callback = wrapcallback(s, disk)
        ),
        storage = isnothing(disk) ? nothing : 0.0
    )
end

"""
    $(TYPEDSIGNATURES)

Dispatch method for integrating a geodesic defined by `s` and `cf`. For configuration of 
this method, see [`BHSetup`](@ref) and [`IntegratorConfig`](@ref).
"""
function integrategeodesic(s::BHSetup, cf::IntegratorConfig; storage = nothing)
    p = GeodesicParams(cf.α, cf.β, s, storage)

    u0 = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    v_temp = (0.0, -1.0, p.θv₀, p.ϕv₀)

    v0 = SVector(null_constrain(u0, v_temp, s.metric), -1.0, p.θv₀, p.ϕv₀)

    integrate(u0, v0, (s.λlow, s.λhigh), p, cf)
end
function integrategeodesic(
    s::BHSetup{CarterBoyerLindquist{T}},
    cf::IntegratorConfig;
    storage = nothing
) where {T}
    p = CarterGeodesicParams(cf.α, cf.β, s, storage)

    u0 = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    integrate(u0, (s.λlow, s.λhigh), p, cf)
end

"""
    $(TYPEDSIGNATURES)
"""
@inline function integrate(x::StaticVector, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{false}(rayintegrator, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(x::StaticVector, v::StaticVector, time_domain, p::GeodesicParams, cf::IntegratorConfig)
    prob = SecondOrderODEProblem{false}(rayintegrator, v, x, time_domain, p)
    solvegeodesic(prob, cf)
end
#= @inline function integrate(x::AbstractVector, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(
    x::AbstractVector,
    time_domain,
    p,
    cf::IntegratorConfig{EnsembleGPUArray}
)
    x = Float32[x...]
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, [changetype(Float32, p)])
    solvegeodesic(prob, cf)
end =#

solvegeodesic(::Any, ::ParallelParams{E,Nothing}) where {E} =
    error("`probfunc` in ParallelParams must be defined.")
function solvegeodesic(prob, cf::ParallelParams{E,P,F}) where {E,P,F}
    solve(
        EnsembleProblem(prob, prob_func = cf.probfunc, safetycopy = false),
        Tsit5(),
        cf.ensemble;
        trajectories = cf.trajectories,
        verbose = cf.verbose,
        maxiters = cf.maxiters,
        abstol = cf.abstol,
        reltol = cf.reltol,
        save_on = cf.save_geodesics,
        callback = cf.callback,
        #dtmax=1.0
    )
end
function solvegeodesic(prob, cf::IntegratorConfig{F}) where {F}
    solve(
        prob,
        Tsit5(),
        ;
        verbose = cf.verbose,
        maxiters = cf.maxiters,
        abstol = cf.abstol,
        reltol = cf.reltol,
        save_on = cf.save_geodesics,
        callback = cf.callback
    )
end


export integrategeodesic, calcgeodesic, rayintegrator, rayintegrator!
