"""
    $(TYPEDSIGNATURES)

In practice, we check ``r > R_0 + \\epsilon``, with ``\\epsilon`` small, and ``r < \\infty``, 
where our effective infinity is hardcoded as `1200.0`.
"""
@inline function chartbounds(r, p)
    on_chart = p.chart_inner_radius < r < 1200.0
    !on_chart
end


"""
    $(TYPEDSIGNATURES)
"""
@inline function is_radial_pot_negative(u, λ, integrator)::Bool
    p = integrator.p
    Vr(p.metric.E, p.L, p.metric.M, p.Q, u[2], p.metric.a) < 0
end


"""
    $(TYPEDSIGNATURES)
"""
@inline function is_angular_pot_negative(u, λ, integrator)::Bool
    p = integrator.p
    Vθ(p.metric.E, p.L, p.Q, p.metric.a, u[3]) < 0
end


"""
    $(TYPEDSIGNATURES)

Updates the [`CarterGeodesicParams`](@ref) stored in the integrator to set the
affine time at which the sign was flipped.
"""
@inline function flip_radial_sign(integrator)
    integrator.p = flip_rsign(integrator.t, integrator.p)
end


"""
    $(TYPEDSIGNATURES)

Updates the [`CarterGeodesicParams`](@ref) stored in the integrator to set the
affine time at which the sign was flipped.
"""
@inline function flip_angular_sign(integrator)
    integrator.p = flip_θsign(integrator.t, integrator.p)
end

"""
    $(TYPEDSIGNATURES)

Creates a callback function which ensures the problem is chart bound. If the disk is not `Nothing`,
additionally checks intersections with the disk.
"""
wrapcallback(s::BHSetup, disk) = DiscreteCallback(
    (u, λ, integrator) -> begin
        # indexing for 2nd order ODE problem has velocity in 1:4, position 5:8
        chartbounds(u[6], integrator.p) || intersect!(integrator, @view(u[5:8]), disk)
    end,
    terminate!
)
wrapcallback(s::BHSetup, disk::Nothing) = DiscreteCallback(
    (u, λ, integrator) -> chartbounds(u[6], integrator.p),
    terminate!
)
function wrapcallback(s::BHSetup{CarterBoyerLindquist{T}}, disk) where {T}
    CallbackSet(
        DiscreteCallback(is_radial_pot_negative, flip_radial_sign),
        DiscreteCallback(is_angular_pot_negative, flip_angular_sign),
        DiscreteCallback(
            (u, λ, integrator) -> chartbounds(u[2], integrator.p) || intersect!(integrator, u, disk),
            terminate!
        )
    )
end
function wrapcallback(s::BHSetup{CarterBoyerLindquist{T}}, disk::Nothing) where {T}
    CallbackSet(
        DiscreteCallback(is_radial_pot_negative, flip_radial_sign),
        DiscreteCallback(is_angular_pot_negative, flip_angular_sign),
        DiscreteCallback(
            (u, λ, integrator) -> chartbounds(u[2], integrator.p),
            terminate!
        )
    )
end
