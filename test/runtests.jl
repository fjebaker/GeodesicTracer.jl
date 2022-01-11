using Test, GeodesicTracer

struct TestMetricParams{T} <: GeodesicTracer.AbstractMetricParams{T}
    test::T
end

GeodesicTracer.constrain(m::TestMetricParams{T}, u, v; μ=μ) where {T} =  sum(v)

SVector = GeodesicTracer.SVector

@testset "Constraints" begin
    tmp = TestMetricParams(1.0)

    @testset "Vectors" begin
        u = Float64[1.0, 1.0, 1.0, 1.0]
        v = Float64[1.0, 2.0, -1.0, -2.0]
        expected = [0.0, 2.0, -1.0, -2.0]

        res_v = GeodesicTracer.constrain_all(tmp, u, v, 0.0)
        @test res_v isa Vector{Float64}
        @test all(res_v .== expected)

        u_static = SVector(u...)
        v_static = SVector(1.0, 2.0, -1.0, -2.0)
        res_v_static = GeodesicTracer.constrain_all(tmp, u_static, v_static, 0.0)
        @test all(res_v_static .== expected)
    end


    @testset "Vector of Vectors" begin 
        us = Vector{Float64}[
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0]
        ]
        vs = Vector{Float64}[
            [1.0, 2.0, -1.0, -2.0],
            [2.0, 2.0, -1.0, -2.0],
            [3.0, 2.0, -1.0, -2.0]
        ]

        res_vs = GeodesicTracer.constrain_all(tmp, us, vs, 0.0)
        for (i, v) in enumerate(res_vs)
            @test v[1] .≈ i - 1.0
        end
    end

end