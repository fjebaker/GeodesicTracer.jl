push!(LOAD_PATH,"src")

using Documenter
using GeodesicTracer

makedocs(
    modules=[GeodesicTracer],
    clean=false,
    sitename="GeodesicTracer Documentation",

    pages = [
        "Home" => "index.md",
        "Coordinates" => "coordinates.md",
        "Integration" => "integration.md",
        "Redshift" => "redshift.md",
        "Rendering" => "rendering.md"
    ]
)

#= 
TODO: once repo is public, so we can use GitHub pages.

deploydocs(
    repo=""
) 
=#