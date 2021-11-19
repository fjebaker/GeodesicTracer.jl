"""
    $(TYPEDEF)
"""
abstract type IntegratorConfig{F,S} end

"""
    $(TYPEDEF)

Used to store the configuration for the solve algorithm when calculating a single geodesic.

# Fields:

$(FIELDS)
"""
@with_kw struct SingleParams{F,S} <: IntegratorConfig{F,S}
    @deftype Float64
    maxiters = 1e4
    reltol = 1e-9
    abstol = 1e-9
    α = 0.0
    β = 0.0

    save_geodesics::Bool = true
    verbose::Bool = false

    callback::F = nothing
    solver::S = Tsit5()
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
@with_kw struct ParallelParams{E,P,F,S} <: IntegratorConfig{F,S}
    @deftype Float64
    maxiters = 1e4
    reltol = 1e-9
    abstol = 1e-9
    α = 0.0
    β = 0.0

    α_end::Float64 = 1.0
    trajectories::Int = 100

    save_geodesics::Bool = true
    verbose::Bool = false

    ensemble::E = nothing
    probfunc::P = nothing
    callback::F = nothing
    solver::S = Tsit5()
end

export SingleParams, ParallelParams
