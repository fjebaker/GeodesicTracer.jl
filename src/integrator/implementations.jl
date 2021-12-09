"""
    $(TYPEDSIGNATURES)
"""
function secondorder_rayintegrator(
    v::StaticVector,
    u::StaticVector,
    p::GeodesicParams{V,T},
    λ
) where {V,T}
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
"""
@inline function integrate(x::StaticVector, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{false}(rayintegrator, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(x, time_domain, p, cf::IntegratorConfig)
    prob = ODEProblem{true}(rayintegrator!, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(x, time_domain, p, cf::IntegratorConfig{EnsembleGPUArray})
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
@inline function integrate(x, v, time_domain, p::GeodesicParams, cf::IntegratorConfig)
    prob = SecondOrderODEProblem{true}(secondorder_rayintegrator!, v, x, time_domain, p)
    solvegeodesic(prob, cf)
end
@inline function integrate(v, x, time_domain, p, cf::IntegratorConfig{EnsembleGPUArray})
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



