# GeodesicTracer documentation

```@meta
CurrentModule = GeodesicTracer
```

## User types

```@docs
BHSetup
```

There is also a special `AbstractVector` type provided.

```@docs
FourVector
```

Note that this is only ever really used for generating function specializations with [`@vec_eq`](@ref), and does not actually get explicitly invoked anymore, since benchmarking tests with `StaticArrays.jl` proved percentage points faster.

It is still a convenient representation of arbitrary four vectors, however, and will be left in for future use.

Some additional points regarding this:
- `FourVector` (and for that matter all of the above structs) pass `isbits` so could be used on the GPU
- `FourVector` may be convenient for creating methods for intersection calculations for arbitrary disk geometries
