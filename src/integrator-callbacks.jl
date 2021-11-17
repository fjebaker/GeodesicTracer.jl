"""
    $(TYPEDSIGNATURES)

In practice, we check ``r > R_0 + \\epsilon``, with ``\\epsilon`` small, and ``r < \\infty``, 
where our effective infinity is hardcoded as `1200.0`.
"""
@inline function chartbounds(u, λ, integrator)
    p = integrator.p

    on_chart = p.chart_inner_radius < u[6] < 1200.0
    !on_chart
end

"""
    $(TYPEDSIGNATURES)

Monitors [`Vr`](@ref) and [`Vθ`](@ref) and flips signs stored repsectively in the 
integrator parameter if they are negative, along with storing the time at which this 
occured.
"""
function signflip(u, λ, integrator)
    p = integrator.p
    metric = p.metric

    if Vr(metric.E, p.L, metric.M, p.Q, u[2], metric.a) < 0
        integrator.p = flip_rsign(λ, p)
    elseif Vθ(metric.E, p.L, p.Q, metric.a, u[3]) < 0
        integrator.p = flip_θsign(λ, p)
    end

    chartbounds(u, λ, integrator)
end

"""
    $(TYPEDSIGNATURES)

Creates a callback function which invoked [`signflip`](@ref) and checks intersections with
the disk given by `disk`.
"""
function wrapcallback(s::BHSetup, disk)
    cb =
        isnothing(disk) ? chartbounds :
        (u, λ, integrator) -> begin
            chartbounds(u, λ, integrator) || intersect!(integrator, u, disk)
        end
    DiscreteCallback(cb, terminate!)
end

function wrapcallback(s::BHSetup{CarterBoyerLindquist{T}}, disk) where {T}
    cb =
        isnothing(disk) ? signflip :
        (u, λ, integrator) -> begin
            signflip(u, λ, integrator) || intersect!(integrator, u, disk)
        end
    DiscreteCallback(cb, terminate!)
end
