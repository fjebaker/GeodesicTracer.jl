abstract type AbstractGeodesicParams{V,T} end

"""
    $(TYPEDSIGNATURES)

Casts [`GeodesicParams`](@ref) of type `{S,T}` to `{S,NewType}`.
"""
function changetype(NewType::Type, p::T) where {T<:AbstractGeodesicParams}
    T(; (f => NewType(getfield(p, f)) for f ∈ fieldnames(T) if f != :storage)...)
end
