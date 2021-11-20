"""
    $(TYPEDEF)

Used to describe the black hole setup for which to trace geodesics and render images. All 
field types are Float64, unless otherwise specified:

$(FIELDS)

The field of view factor relates as to how the image pixels at ``x`` and ``y`` 
from the center of the image are mapped into impact parameters ``\\alpha`` and ``\\beta``, 
via

```math
\\left( \\alpha, \\beta \\right)
    = \\frac{1}{f} \\left(x, y\\right),
```

where ``f`` is the field of view factor.
"""
@with_kw mutable struct BHSetup{M}
    @deftype Float64
    "Metric structure."
    metric::M = CarterBoyerLindquist()

    "Azimuthal angle of observer."
    ϕ₀ = 0.0             # azimuthal angle
    "Inclination of observer."
    θ₀ = deg2rad(90.0)   # inclination
    "Initial position of observer."
    r₀ = 1000.0

    "Lower bounds on the integrator time (there is no practical reason for this not to be 0.0."
    λlow = 0.0
    "Upper bounds on the integrator time."
    λhigh = 2000.0

    "Width of the rendered image in pixels (`Int`)."
    img_width::Int = 300
    "Height of the rendered image in pixels (`Int`)."
    img_height::Int = 180
    "Field of view factor."
    fov_factor = 3.0    # field of view factor
end

export BHSetup
