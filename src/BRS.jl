module GeodesicTracer

using DifferentialEquations
using StaticArrays
using DiffEqGPU

import Base: getindex, setindex!, size, length

using DocStringExtensions
using Parameters
using Accessors

include("bh-setup.jl")
include("geodesic-params.jl")
include("integrator-config.jl")

include("fanton.jl")

include("fourvector.jl")
include("utils.jl")
include("disks.jl")

include("coordinates.jl")
include("redshift.jl")
include("integrator.jl")
include("value-functions.jl")
include("render.jl")

end
