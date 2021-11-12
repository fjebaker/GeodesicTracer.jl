"""
    $(TYPEDEF)
"""
abstract type IntegratorConfig{F} end

"""
    $(TYPEDEF)

Used to store the configuration for the solve algorithm when calculating a single geodesic.

# Fields:

$(FIELDS)
"""
@with_kw struct SingleParams{F} <: IntegratorConfig{F}
    @deftype Float64
    maxiters = 1e4
    reltol = 1e-8
    abstol = 1e-8
    α = 0.0
    β = 0.0

    save_geodesics::Bool = true
    verbose::Bool = false

    callback::F = nothing
end


"""
    $(TYPEDEF)

Used to store the configuration for the solve algorithm when calculating a multiple 
geodesics. The `E` parameter parameter is the ensemble algorithm to use from
[DifferentialEquations.jl](https://diffeq.sciml.ai/stable/features/ensemble/).

Commonly used:

  - EnsembleThreads
  - EnsembleGPUArray (from DiffEqGPU.jl)

# Fields:

$(FIELDS)
"""
@with_kw struct ParallelParams{E,P,F} <: IntegratorConfig{F}
    @deftype Float64
    maxiters = 1e4
    reltol = 1e-8
    abstol = 1e-8
    α = 0.0
    β = 0.0

    α_end::Float64 = 1.0
    trajectories::Int = 100

    save_geodesics::Bool = true
    verbose::Bool = false

    ensemble::E = nothing
    probfunc::P = nothing
    callback::F = nothing
end

export SingleParams, ParallelParams
