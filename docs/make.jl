push!(LOAD_PATH,"src")

using Documenter
using GeodesicTracer

makedocs(
    modules=[GeodesicTracer],
    clean=false,
    sitename="GeodesicTracer Documentation",

    pages = [
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Setup & Configuration" => "bh_setup.md",
        "Integration" => [
            "1st Order" => "first_order_integration.md",
            "2nd Order" => "second_order_integration.md"
        ],
        "Rendering" => "rendering.md",
        "Value Functions" => [
            "Generic" => "value_functions.md",
            "Redshift" => "redshift.md"
        ],
    ]
)

#= 
TODO: once repo is public, so we can use GitHub pages.

deploydocs(
    repo=""
) 
=#