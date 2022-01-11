# single position and single velocity
function __tracegeodesics(
    m::AbstractMetricParams{T}, 
    init_pos::AbstractVector{T}, 
    init_vel::AbstractVector{T},
    time_domain::Tuple{T,T},
    solver
    ; 
    μ,
    kwargs...
    ) where {T}
    prob = integrator_problem(m, init_pos, init_vel, time_domain)
    integrate_prob(m, solver, prob; kwargs...)
end

# columnar of positions and velocity
function __tracegeodesics(
    m::AbstractMetricParams{T}, 
    init_positions::AbstractArray{T, 2}, 
    init_velocities::AbstractArray{T, 2},
    time_domain::Tuple{T,T},
    solver
    ;
    kwargs...
    ) where {T}
    prob = integrator_problem(m, @view(init_positions[:, 1]), @view(init_velocities[:, 1]), time_domain)
    ens_prob = EnsembleProblem(
        prob,
        prob_func = (prob, i, repeat) -> integrator_problem(m, @view(init_positions[:, i]), @view(init_velocities[:, i]), time_domain),
        safetycopy = false
    )
    integrate_prob(m, solver, ens_prob; ensemble = EnsembleThreads(), trajectories=length(init_positions), kwargs...)
end

function integrate_prob(m::AbstractMetricParams{T}, solver, prob; callbacks, kwargs...) where {T}
    cbs = create_callback_set(m, callbacks)
    solve_prob_with_cbs(solver, prob, cbs; kwargs...)
end

function solve_prob_with_cbs(solver, prob, cbs; ensemble, solver_opts...)
    solve(
        prob,
        solver,
        ensemble
        ;
        callback=cbs,
        solver_opts...
    )
end
