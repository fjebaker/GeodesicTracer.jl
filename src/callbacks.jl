function metric_callback(m::AbstractMetricParams{T}) where {T}
    min_r = inner_radius(m)
    DiscreteCallback(
        (u, λ, integrator) -> u[6] ≤ min_r * 1.1 || u[6] > 1200.0,
        terminate!
    )
end

function create_callback_set(m::AbstractMetricParams{T}, cb::Nothing) where {T}
    metric_callback(m)
end

function create_callback_set(m::AbstractMetricParams{T}, cb::Tuple) where {T}
    CallbackSet(
        metric_callback(m),
        cb...
    )
end