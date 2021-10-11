using Test, GeodesicTracer

@testset "Types and Constructors" begin
    s = BHSetup()
    p = GeodesicTracer.GeodesicParams(0.0, 0.0, s)

    @test typeof(GeodesicTracer.changetype(Float32, p).a) == Float32
end

@testset "Coordinate Consistency" begin
    s = BHSetup()

    p = GeodesicTracer.GeodesicParams(4.0, 0.0, s)
    @test GeodesicTracer.rms(s) == p.rms
    @test GeodesicTracer.rms(s) == GeodesicTracer.rms(p)

    s.a = -1.0
    p = GeodesicTracer.GeodesicParams(4.0, 0.0, s)
    @test GeodesicTracer.rms(s) == p.rms
    @test GeodesicTracer.rms(s) == GeodesicTracer.rms(p)

    s.a = 1.0
    p = GeodesicTracer.GeodesicParams(4.0, 0.0, s)
    @test GeodesicTracer.rms(s) == p.rms
    @test GeodesicTracer.rms(s) == GeodesicTracer.rms(p)

    p2 = GeodesicTracer.flip_rsign(9.0, p)
    @test p2.r_sign == - p.r_sign
    @test p2.λr_change == 9.0

    p2 = GeodesicTracer.flip_θsign(9.0, p)
    @test p2.θ_sign == - p.θ_sign
    @test p2.λθ_change == 9.0
end 

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