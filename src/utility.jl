"""
    map_to_velocity(m::AbstractMetricParams{T}, u, α, β)

Map impact parameters `α` and `β` to a velocity vector at some position `u` in the given metric `m`.
"""
function map_to_velocity(m::AbstractMetricParams{T}, u::AbstractVector{T}, α, β) where {T}
    [0.0, alpha_beta_to_vel(m, u, α, β)...]
end

function map_to_velocity(m::AbstractMetricParams{T}, u::SVector{S,T}, α, β) where {S,T}
    SVector{S,T}(0.0, alpha_beta_to_vel(m, u, α, β)...)
end

function map_to_velocity(m::AbstractMetricParams{T}, u, α::AbstractVector{T}, β::AbstractVector{T}) where {S,T}
    curried(_u, _α, _β) = map_to_velocity(m, _u, _α, _β)
    curried.(u, α, β)
end

function alpha_beta_to_vel(m::AbstractMetricParams{T}, u, α, β) where {T}
    reg = u[2]^2
    -1, α / reg, β / reg
end

export map_to_velocity