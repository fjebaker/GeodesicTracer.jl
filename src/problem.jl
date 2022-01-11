function integrator_problem(
    m::AbstractMetricParams{T},
    pos::StaticVector{T}, vel::StaticVector{T}, time_domain
    ) where {T}
    SecondOrderODEProblem{false}(vel, pos, time_domain, m) do u, v, p, λ
        SVector(geodesic_eq(u, v, p)...)
    end
end

function integrator_problem(
    m::AbstractMetricParams{T},
    pos::AbstractVector{T}, vel::AbstractVector{T}, time_domain
    ) where {T}
    SecondOrderODEProblem{true}(vel, pos, time_domain, m) do dv, u, v, p, λ
        dv .= geodesic_eq(u, v, p)
    end
end