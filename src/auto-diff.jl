# for use with static, axis-symmetric metrics
# create new abstract type for easy re-definition
abstract type AbstractAutoDiffMetricParams{T} <: AbstractMetricParams{T} end
abstract type AbstractAutoDiffStaticAxisSymmetricParams{T} <: AbstractAutoDiffMetricParams{T} end

# interface
metric_components(m::AbstractAutoDiffMetricParams{T}, r, θ) where {T} = error("Not implemented for $(typeof(m)).")

# this function could just be used to generate the metric matrix, since
# the jacobian for static, axis-symmetric metrics has the same non-zero components
# of the christoffel connections, so we can bypass having to use Tullio.jl and therefore
# having to turn the jacobian into a matrix
function to_matrix(g_comp, g_jac_components)
    g = @SMatrix [
        g_comp[1] 0.0 0.0 g_comp[5] ;
        0.0 g_comp[2] 0.0 0.0 ;
        0.0 0.0 g_comp[3] 0.0 ;
        g_comp[5] 0.0 0.0 g_comp[4]
    ]

    j1 = @SMatrix [
        0.0 0.0 0.0 0.0 ;
        0.0 0.0 0.0 0.0 ;
        0.0 0.0 0.0 0.0 ;
        0.0 0.0 0.0 0.0
    ]

    dr = g_jac_components[1]
    dθ = g_jac_components[2]

    j2 = @SMatrix [
        dr[1] 0.0 0.0 dr[5] ;
        0.0 dr[2] 0.0 0.0 ;
        0.0 0.0 dr[3] 0.0 ;
        dr[5] 0.0 0.0 dr[4]
    ]

    j3 = @SMatrix [
        dθ[1] 0.0 0.0 dθ[5] ;
        0.0 dθ[2] 0.0 0.0 ;
        0.0 0.0 dθ[3] 0.0 ;
        dθ[5] 0.0 0.0 dθ[4]
    ]

    return g, (j1, j2, j3, j1)
end

function christoffel_components(g_jac, ginv)
    @tullio Γ[i, k, l] :=
        1 / 2 *
        ginv[i, m] *
        # l is the index of the derivative
        (g_jac[l][m, k] + g_jac[k][m, l] - g_jac[m][k, l])
end

calc_geodesic_equation(Γ, v) = @tullio δxδλ[i] := -v[j] * v[k] * Γ[i, j, k]

# should eventually look into using DiffResults.jl to do both the jacobian and metric 
# evaluation in a single go, but this works for now
# the difficulty is that the output of the metric params is of greater dimension than
# the input, so the straight forward implementation isn't quite so good
# also should check if ForwardDiff.jl can work directly with static arrays
@inbounds function geodesic_eq(m::AbstractAutoDiffMetricParams{T}, u, v) where {T}
    g_comp = metric_components(m, u[2], u[3])
    jacs = ForwardDiff.jacobian(
        x -> metric_components(m, x[1], x[2]), 
        @view(u[2:3])
    )      
    j1 = SVector{5, Float64}(jacs[:, 1])
    j2 = SVector{5, Float64}(jacs[:, 2])
    
    g, jac = to_matrix(g_comp, (j1, j2))
    ginv = inv(g)
    Γ = christoffel_components(jac, ginv)

    SVector{4, Float64}(calc_geodesic_equation(Γ, v))
end

export AbstractAutoDiffStaticAxisSymmetricParams