"""
    $(TYPEDSIGNATURES)

Creates specialisations for the function `f` given in `expr`. 

Wrapping

```julia
@vec_eq f(r, a, θ) = ...
```

generates definitions for

```julia
f(r, a, θ) = ...
f(x::FourVector, p::GeodesicParams) = f(x.r, p,a, x.θ)
f(x::AbstractVector, p::GeodesicParams) = f(x[2], p.a, x[3])
```

"""
macro vec_eq(expr)
    # get function name
    f = expr.args[1]
    name = f.args[1]

    # store arguments not present in FourVector or GeodesicParams
    new_args::Array{Union{Expr,Symbol}} = [:(p::GeodesicParams)]

    # make appropriate argument substitution
    arguments = []
    for arg in f.args[2:end]
        if hasfield(FourVector, arg)
            arg = :(x.$arg)
        elseif hasfield(GeodesicParams, arg)
            arg = :(p.$arg)
        else
            push!(new_args, arg)
        end
        push!(arguments, arg)
    end

    # map the symbolic expression to a number index version 
    s_arguments = map(
        i -> begin
            if i isa Expr && i.args[1] == :x
                index = findfirst(==(i.args[2].value), fieldnames(FourVector))
                i = :(x[$index])
            end
            i
        end,
        arguments
    )

    # build function definitions
    quote
        Base.@__doc__ @inline $(expr)
        @inline @inbounds $(name)(x::FourVector, $(new_args...)) =
            $(name)($(arguments...))
        @inline @inbounds $(name)(x::AbstractVector, $(new_args...)) =
            $(name)($(s_arguments...))
    end |> esc
end
