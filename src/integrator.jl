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
function rayintegrator(u, p, λ)
    δ(u, p)
end

"""
    $(TYPEDSIGNATURES)

Inplace variant of [`rayintegrator`](@ref).
"""
function rayintegrator!(du, u, p, λ)
    δ!(du, u, p)
end
function rayintegrator!(du, u, p::AbstractArray{GeodesicParams}, λ)
    # gpu version requires parameter unpacking
    δ!(du, u, p[1])
end

"""
    $(TYPEDSIGNATURES)

Monitors [`Vr`](@ref) and [`Vθ`](@ref) and flips signs stored repsectively in the 
integrator parameter if they are negative, along with storing the time at which this 
    occured.

Furthermore, is used in the `DiscreteCallback` of the integrator to check if the geodesic is 
still in the local chart.

In practice, we check ``r > R_0 + \\epsilon``, with ``\\epsilon`` small, and ``r < \\infty``, 
where our effective infinity is hardcoded as `1200.0`.
"""
function signflip(u, λ, integrator)
    p = integrator.p

    if Vr(u, p) < 0
        integrator.p = flip_rsign(λ, p)
    elseif Vθ(u, p) < 0
        integrator.p = flip_θsign(λ, p)
    end

    on_chart = p.R₀ < u[2] < 1200.0
    !on_chart
end

"""
    $(TYPEDSIGNATURES)

Creates a callback function which invoked [`signflip`](@ref) and checks intersections with
the disk given by `disk`.
"""
function wrapcallback(disk)
    (u, λ, integrator) -> begin
        signflip(u, λ, integrator) || intersect!(integrator, u, disk)
    end
end


"""
    $(TYPEDSIGNATURES)
"""
function newparams(p::GeodesicParams, θ, r, α, β, δα)::GeodesicParams
    l, q = LQ(p.M, r, p.a, θ, α + δα, β)
    @set(@set(p.L = l).Q = q)
end


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
            callback = callback = DiscreteCallback(
                isnothing(disk) ? signflip : wrapcallback(disk),
                terminate!
            )
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
    δα = (α_range[2] - α_range[1]) / num
    probfunc(prob, i, repeat) =
        remake(prob, p = newparams(prob.p, prob.u0[3], prob.u0[2], α_range[1], β, i * δα))
    integrategeodesic(
        s,
        ParallelParams(
            trajectories = num,
            β = β,
            α = α_range[1],
            α_end = α_range[2],
            save_geodesics = save_geodesics,
            ensemble = EnsembleThreads(),
            probfunc = probfunc,
            callback = DiscreteCallback(
                isnothing(disk) ? signflip : wrapcallback(disk),
                terminate!
            )
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
    integrate(s, p, cf)
end

# fallback
@inline function integrate(s::BHSetup, p, cf::IntegratorConfig)
    x = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    integrate(x, (s.λlow, s.λhigh), p, cf)
end

# GPU
@inline function integrate(s::BHSetup, p, cf::ParallelParams{EnsembleGPUArray})
    # make everything float 32 for GPU
    x = Float32[0.0f0, s.r₀, s.θ₀, s.ϕ₀]
    integrate(x, (Float32(s.λlow), Float32(s.λhigh)), [changetype(Float32, p)], cf)
end

@inline function integrate(x::StaticVector, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{false}(rayintegrator, x, time_domain, p)
    #(prob, cf)
    solvegeodesic(prob, cf)
end
@inline function integrate(x::AbstractVector, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, p)
    solvegeodesic(prob, cf)
end

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
        callback = cf.callback
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
