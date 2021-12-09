"""
    $(TYPEDSIGNATURES)
"""
@inline function render!(row, α_range, β, d, s, valuefunc::ValueFunction; solver=Tsit5())
    simsols = calcgeodesic(α_range, s.img_width, β, s, disk = d, save_geodesics = false, solver=solver)
    Threads.@threads for i = 1:s.img_width
        row[i] = valuefunc(simsols[i], s, d)
    end
end

"""
    $(TYPEDSIGNATURES)

Renders an image (currently only a redshift map) of an `AccretionDisk` a black hole 
described by `s`, a [`BHSetup`](@ref).

This is achieved by tracing a single geodesic for each pixel in the image, and calculating 
the redshift if it intercepts the disk.

Allocates and returns a 2-dimensional `Float64` array.
"""
function render(d, s; valuefunc = geometry, kwargs...)
    # preallocate image output
    image = zeros(Float64, (s.img_height, s.img_width))

    render!(image, d, s, valuefunc; kwargs...)
    image
end
render(s; kwargs...) = render(nothing, s; kwargs...)


"""
    $(TYPEDSIGNATURES)

Inplace version of [`renderdisk`](@ref), where the image is allocated by the user.
"""
function render!(image, d, s, valuefunc::ValueFunction; kwargs...)
    y_mid = s.img_height ÷ 2
    x_mid = s.img_width ÷ 2

    @inbounds for Y = 1:s.img_height
        # generate space of α to integrate
        # have to use a slight 0.01 offset to avoid integrating α=0.0 geodesics
        α_range = ((0.01 - x_mid) / s.fov_factor, (0.01 + x_mid) / s.fov_factor)
        render!(@view(image[Y, :]), α_range, (Y - y_mid) / s.fov_factor, d, s, valuefunc; kwargs...)
    end
end

export render, render!
