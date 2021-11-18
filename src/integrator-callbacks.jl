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

Monitors [`Vr`](@ref) and [`Vθ`](@ref) and flips signs stored repsectively in the 
integrator parameter if they are negative, along with storing the time at which this 
occured.
"""
@inline function signflip(u, λ, integrator)
    p = integrator.p
    metric = p.metric

    if Vr(metric.E, p.L, metric.M, p.Q, u[2], metric.a) < 0
        integrator.p = flip_rsign(λ, p)
    elseif Vθ(metric.E, p.L, p.Q, metric.a, u[3]) < 0
        integrator.p = flip_θsign(λ, p)
    end
end

@inline function is_radial_pot_negative(u, λ, integrator)
    p = integrator.p
    Vr(p.metric.E, p.L, p.metric.M, p.Q, u[2], p.metric.a) < 0
end

@inline function is_angular_pot_negative(u, λ, integrator)
    p = integrator.p
    Vθ(p.metric.E, p.L, p.Q, p.metric.a, u[3]) < 0
end

@inline function flip_radial_sign(integrator)
    integrator.p = flip_rsign(integrator.t, integrator.p)
end

@inline function flip_angular_sign(integrator)
    integrator.p = flip_θsign(integrator.t, integrator.p)
end


"""
    $(TYPEDSIGNATURES)

Creates a callback function which invoked [`signflip`](@ref) and checks intersections with
the disk given by `disk`.
"""
function wrapcallback(s::BHSetup, disk)
    cb =
        isnothing(disk) ? (u, λ, integrator) -> chartbounds(u[6], integrator.p) :
        (u, λ, integrator) -> begin
            # indexing for 2nd order ODE problem has velocity in 1:4, position 5:8
            chartbounds(u[6], integrator.p) || intersect!(integrator, @view(u[5:8]), disk)
        end
    #CallbackSet(
    DiscreteCallback(cb, terminate!)
    #DiscreteCallback(
    #    (u, λ, integrator) -> u[6] < 5 * integrator.p.chart_inner_radius,
    #    adjust_time_step
    #    )
    #)
end

function wrapcallback(s::BHSetup{CarterBoyerLindquist{T}}, disk) where {T}
    chart_callback =
        isnothing(disk) ? (u, λ, integrator) -> begin 
            chartbounds(u[2], integrator.p)
        end :
        (u, λ, integrator) -> begin
            chartbounds(u[2], integrator.p) || intersect!(integrator, u, disk)
        end
    CallbackSet(
        DiscreteCallback(is_radial_pot_negative, flip_radial_sign),
        DiscreteCallback(is_angular_pot_negative, flip_angular_sign),
        DiscreteCallback(chart_callback, terminate!)
    )
end
