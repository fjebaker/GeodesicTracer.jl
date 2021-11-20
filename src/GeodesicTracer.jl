module GeodesicTracer

using DifferentialEquations
using StaticArrays
using DiffEqGPU

import Base: getindex, setindex!, size, length

using DocStringExtensions
using Parameters
using Accessors

using ComputedGeodesicEquations

# re-export metric structures
export BoyerLindquist, EddingtonFinkelstein
# include additional metric structures
include("carter-method/carter-boyer-lindquist.jl")

# setup functions
include("bh-setup.jl")

# geodesic parameters
include("abstract-geodesic-params.jl")
include("geodesic-eq-method/geodesic-params.jl")
include("carter-method/carter-params.jl")

# integrator setup
include("integrator-config.jl")
include("integrator-callbacks.jl")
include("integrator.jl")

# coordinate functions
include("coordinates.jl")
include("carter-method/carter-quantities.jl")
include("carter-method/fanton.jl")
include("redshift.jl")

# disks and rendering
include("disks.jl")
include("value-functions.jl")
include("render.jl")

end
