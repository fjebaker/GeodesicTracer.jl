abstract type AccretionDisk end

"""
    $(TYPEDEF)

Geometrically thin, optically thick accretion disk.

$(FIELDS)
"""
@with_kw struct GeometricDisk{T} <: AccretionDisk
    "Inner radius of the disk"
    r_inner::T = 10.0
    "Outer radius of the disk"
    r_outer::T = 50.0
    "Inclination of the disk towards the observer"
    inclination::T = π / 2
end

"""
    $(TYPEDSIGNATURES)

Calculates whether the point ``u`` intersects the disk, projected from the previous point
(state stored in `integrator.p`, a [`GeodesicParams`](@ref)).

Returns `true` if the point intersects, `false` otherwise. Updates the `out` field of `p`.

# Details

The disk is assumed to be the disk lying on a plane given by the normal vector ``\\vec{n} = 
(\\cos \\alpha, 0, \\sin \\alpha)``. This
vector describes a plane tilted along the ``y`` axis towards the positive ``x`` axis by 
some angle ``\\alpha``.

A line connecting the points ``\\vec{p}_1`` and ``\\vec{p}_2`` will intercept a plane if

```math
\\left( \\vec{n} \\cdot \\vec{p}_1 \\right) 
    \\left( \\vec{n} \\cdot \\vec{p}_2 \\right) \\leq 0,
```

i.e. if the sign of the two dot products differs.

Using the coordinate transformations ``x = r \\sin \\theta \\cos \\phi`` and ``z = r \\cos 
\\phi``, along with the above defined ``\\vec{n}``, this condition may be imposed for 
spherical coordinates:

```math
r_1 r_2 
\\left( 
    \\cos \\alpha \\sin \\theta_1 \\cos \\phi_1 + \\sin \\alpha \\cos \\theta_1
\\right)
\\left( 
    \\cos \\alpha \\sin \\theta_2 \\cos \\phi_2 + \\sin \\alpha \\cos \\theta_2
\\right)
    \\leq 0.
```
"""
function intersect!(integrator, u, d::GeometricDisk)
    sinα = sin(d.inclination)
    cosα = cos(d.inclination)
    p = integrator.p
    r = u[2]
    θ = u[3]
    ϕ = u[4]

    s = r * (cosα * sin(θ) * cos(ϕ) + sinα * cos(θ))
    if (s * p.storage) ≤ 0 && d.r_inner ≤ r ≤ d.r_outer
        return true
    end
    # update the stored value
    integrator.p = @set(p.storage = s)

    false
end

export GeometricDisk
