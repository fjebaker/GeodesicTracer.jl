module GeodesicTracer

using LinearAlgebra
import Base: getindex, setindex!, size, length

using DifferentialEquations
using StaticArrays
using DiffEqGPU
using RecursiveArrayTools

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
include("geodesic-eq-method/geodesic-params.jl")
include("carter-method/carter-params.jl")

# integrator setup
include("integrator/configurations.jl")
include("integrator/callbacks.jl")
include("integrator/implementations.jl")
include("integrator/interface.jl")

# coordinate functions
include("coordinates.jl")
include("carter-method/carter-quantities.jl")
include("carter-method/fanton.jl")

# disks and rendering
include("disks.jl")
include("value-functions.jl")
include("render.jl")

include("redshift.jl")

end
