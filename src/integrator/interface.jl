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
function integrategeodesic(s::BHSetup{M}, cf::IntegratorConfig; storage = nothing) where {M}
    p = GeodesicParams(cf.α, cf.β, s, storage)

    x = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    v_temp = (0.0, -1.0, p.θv₀, p.ϕv₀)

    v = SVector(null_constrain(x, v_temp, s.metric), -1.0, p.θv₀, p.ϕv₀)

    integrate(v, x, (s.λlow, s.λhigh), p, cf)
end
"""
Specialisation for 1st order case.
"""
function integrategeodesic(
    s::BHSetup{CarterBoyerLindquist{T}},
    cf::IntegratorConfig;
    storage = nothing
) where {T}
    p = CarterGeodesicParams(cf.α, cf.β, s, storage)

    u0 = SVector(0.0, s.r₀, s.θ₀, s.ϕ₀)
    integrate(u0, (s.λlow, s.λhigh), p, cf)
end

export integrategeodesic, calcgeodesic