const __four_vec_keys = Dict{Int,Symbol}(1 => :t, 2 => :r, 3 => :θ, 4 => :ϕ)

"""
    $(TYPEDEF)

$(FIELDS)

A radial four vector implementation, providing the `AbstractVector` interface.
"""
struct FourVector{T} <: AbstractVector{T}
    t::T
    r::T
    θ::T
    ϕ::T

    function FourVector(v::AbstractVector{T}) where {T}
        new{T}(v...)
    end
    function FourVector(v::AbstractMatrix{T}) where {T}
        new{T}(v...)
    end
    function FourVector(t::NTuple{4,T}) where {T}
        new{T}(t...)
    end
    function FourVector(t::T, r::T, θ::T, ϕ::T) where {T}
        new{T}(t, r, θ, ϕ)
    end
    function FourVector{T}(i::T...) where {T}
        FourVector(i...)
    end
end

Base.getindex(x::FourVector, i::Int) = getfield(x, __four_vec_keys[i])
Base.setindex!(x::FourVector{T}, value::T, i::Int) where {T} =
    setfield!(x, __four_vec_keys[i], value)
Base.size(::FourVector) = (4,)
Base.length(::FourVector) = 4

export FourVector
