"""
    $(TYPEDSIGNATURES)
"""
function secondorder_rayintegrator(v::StaticVector, u::StaticVector, p::GeodesicParams{V,T}, λ) where {V,T}
    SVector(geodesic_eq(u, v, p.metric)...)
end


"""
    $(TYPEDSIGNATURES)

Inplace variant of [`secondorder_rayintegrator`](@ref).
"""
function secondorder_rayintegrator!(dv, v, u, p::GeodesicParams{V,T}, λ) where {V,T}
    dv .= geodesic_eq(u, v, p.metric)
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
    s::BHSetup{T};
    save_geodesics = true,
    disk = nothing,
    solver = Tsit5()
) where {T}
    # single geodesic method
    integrategeodesic(
        s,
        SingleParams(
            α = α,
            β = β,
            save_geodesics = save_geodesics,
            callback = wrapcallback(s, disk),
            solver = solver
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
    s::BHSetup{T};
    save_geodesics = true,
    disk = nothing,
    solver = Tsit5(),
    ensemble = EnsembleThreads()
) where {T}

    integrategeodesic(
        s,
        ParallelParams(
            trajectories = num,
            β = β,
            α = α_range[1],
            α_end = α_range[2],
            save_geodesics = save_geodesics,
            ensemble = ensemble,
            probfunc = makeprobfunc(s, α_range, β, num),
            callback = wrapcallback(s, disk),
            solver = solver
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

    x = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    v_temp = (0.0, -1.0, p.θv₀, p.ϕv₀)

    v = SVector(null_constrain(x, v_temp, s.metric), -1.0, p.θv₀, p.ϕv₀)

    integrate(v, x, (s.λlow, s.λhigh), p, cf)
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
@inline function integrate(x, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(
    x,
    time_domain,
    p,
    cf::IntegratorConfig{EnsembleGPUArray}
)
    x = Float32[x...]
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, [changetype(Float32, p)])
    solvegeodesic(prob, cf)
end
# second order variants -- TODO: use metaprogramming to generate these?
@inline function integrate(
    x::StaticVector,
    v::StaticVector,
    time_domain,
    p::GeodesicParams,
    cf::IntegratorConfig
)
    prob = SecondOrderODEProblem{false}(secondorder_rayintegrator, v, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(
    x,
    v,
    time_domain,
    p::GeodesicParams,
    cf::IntegratorConfig
)
    prob = SecondOrderODEProblem{true}(secondorder_rayintegrator!, v, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(
    v,
    x,
    time_domain,
    p,
    cf::IntegratorConfig{EnsembleGPUArray}
)
    v_float32 = Float32[v...]
    x_float32 = Float32[x...]
    prob = SecondOrderODEProblem{true}(
        secondorder_rayintegrator!,
        v_float32,
        x_float32,
        time_domain,
        [changetype(Float32, p)]
    )
    solvegeodesic(prob, cf)
end

solvegeodesic(::Any, ::ParallelParams{E,Nothing}) where {E} =
    error("`probfunc` in ParallelParams must be defined.")
function solvegeodesic(prob, cf::ParallelParams{E,P,F,S}) where {E,P,F,S}
    solve(
        EnsembleProblem(prob, prob_func = cf.probfunc, safetycopy = false),
        cf.solver,
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
function solvegeodesic(prob, cf::IntegratorConfig{F,S}) where {F,S}
    solve(
        prob,
        cf.solver,
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
