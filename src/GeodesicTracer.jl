module GeodesicTracer

using SciMLBase
using OrdinaryDiffEq
using DiffEqCallbacks

using StaticArrays
using DocStringExtensions
using Parameters

import GeodesicBase: AbstractMetricParams, geodesic_eq, constrain, on_chart, inner_radius

import ForwardDiff
import Tullio: @tullio

include("callbacks.jl")
include("problem.jl")
include("tracer.jl")
include("constraints.jl")
include("utility.jl")
include("auto-diff.jl")

"""
    tracegeodesics(
        m::AbstractMetricParams{T}, 
        init_positions, init_velocities, 
        time_domain::Tuple{T,T}
        ; 
        μ = 0.0f0, 
        callbacks=Nothing,
        solver=Tsit5(),
        solver_opts...
    ) 

Integrate a geodesic for metric parameterised by `m`, for some initial positions and velocities.
The positions and velocities may be

  - a single position and velocity in the form of a vector of numbers
  - a collection of positions and velocities, as either a vector of vectors, or as a matrix

The matrix specification reads each corresponding column as the initial position and velocity. When a collection of
positions and velocities is supplied, this method dispatched `EnsembleProblem`, offering `ensemble` as a `solver_opts`,
specifying the ensemble method to use.

`solver_opts` are the common solver options in DifferentialEquations.
"""
function tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions,
    init_velocities,
    time_domain::Tuple{T,T};
    solver = Tsit5(),
    μ = 0.0,
    solver_opts...
) where {T}
    __tracegeodesics(
        m,
        init_positions,
        # ensure everything already normalised
        constrain_all(m, init_positions, init_velocities, T(μ)),
        time_domain,
        solver;
        callback = nothing,
        abstol = 1e-8,
        reltol = 1e-8,
        solver_opts...
    )
end

export tracegeodesics

end # module
