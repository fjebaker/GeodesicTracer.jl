using Test, GeodesicTracer

@testset "Geodesic Integrations" begin
    s = BHSetup()

    @testset "Single Geodesic Integrations" begin
        sol = calcgeodesic(4.0, 0.1, s)
        @test sol.retcode == :MaxIters
    end

    @testset "Ensemble Geodesic Integrations" begin
        simsol = calcgeodesic((0.1, 4.0), 10, 0.1, s)
        for sol in simsol
            @test sol.retcode == :MaxIters
        end
    end

    d = GeometricDisk(r_inner=GeodesicTracer.rms(s), inclination=deg2rad(80.0))
    @testset "Disk Collision Integrations" begin
        sol = calcgeodesic(4.0, 0.1, s, disk=d)
        @test sol.retcode == :MaxIters
        
        simsol = calcgeodesic((0.1, 4.0), 10, 0.1, s, disk=d)
        for sol in simsol
            @test sol.retcode == :MaxIters
        end
    end

end