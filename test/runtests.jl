using SoilMoisture
using Test

@testset "SoilMoisture.jl" begin
    # Parameters
    s1 = 0.10
    s2 = 0.20
    s3 = 0.45
    s4 = 1.00
    sh = 0.12
    sw = 0.25
    sstar = 0.50
    sfc = 0.65
    emax = 0.50
    ew = emax / 10
    b = 4.0
    ks = 20
    n = 0.50
    zr = 40.0
    dt = 1 / 48
    rain1 = 0.0
    rain2 = 10.0
    rain3 = zeros(48 * 10)
    params = (
        sh = sh,
        sw = sw,
        sstar = sstar,
        sfc = sfc,
        emax = emax,
        ew = ew,
        b = b,
        ks = ks,
        n = n,
        zr = zr,
        dt = dt
             )

    # Tests
    @test rainfall_poisson(1, 0.0, 0.0)[1] == 0.0
    @test rainfall_poisson(1, 1.0, 1.0)[1] > 0.0
    @test evapotranspiration(s1, sh, sw, sstar, emax, ew) == 0.0 
    @test evapotranspiration(s4, sh, sw, sstar, emax, ew) == emax 
    @test round(evapotranspiration(s2, sh, sw, sstar, emax, ew), digits=2) == 0.03 
    @test round(evapotranspiration(s3, sh, sw, sstar, emax, ew), digits=2) == 0.41 
    @test leakage(s1, sfc, b, ks) == 0.0
    @test leakage(s4, sfc, b, ks) == ks
    @test water_loss(s1, sh, sw, sstar, sfc, emax, ew, b, ks) == 0.0
    @test water_loss(s4, sh, sw, sstar, sfc, emax, ew, b, ks) == (ks + emax)
    @test soil_water_balance(rain1, s1, sh, sw, sstar, sfc, b, ks, n, zr, emax, ew, dt)[1] == s1
    @test soil_water_balance(rain2, s4, sh, sw, sstar, sfc, b, ks, n, zr, emax, ew, dt)[4] == rain2
    res1 = solve_swb(rain3, s3, params)
    res2 = dt2daily(res1)
    @test sum(res1.Q) == 0.0
    @test size(res2)[1] == 11

end
