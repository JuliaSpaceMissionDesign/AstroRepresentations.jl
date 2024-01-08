
@testset "Angle conversions" begin 

    @test_throws ArgumentError("eccentricity must be greater than 0") man_to_tra(0.0, -1.0)

    # ---
    # elliptic
    # Vallado, D. A. (2001). Fundamentals of astrodynamics and applications, pg. 233
    M = 235.4 * pi/180  # rad 
    e = 0.4 

    expected_E = 3.848_661_46 # rad
    expected_nu = eca_to_tra(expected_E, e)

    # get true anomaly
    nu = man_to_tra(M, e)
    @test abs(expected_nu - nu)*180/pi < 1e-2

    # going back to mean anomaly
    M_ = tra_to_man(nu, e)
    @test abs(mod(M - M_, 2*pi))*180/pi < 1e-8

    # --- 
    # hyperbolic 
    # Vallado, D. A. (2001). Fundamentals of astrodynamics and applications, pg. 239
    M = 235.4 * pi/180  # rad 
    e = 2.4 

    expected_H = 1.601_396_28 # rad
    expected_nu = hya_to_tra(expected_H, e)

    # get true anomaly
    nu = man_to_tra(M, e)
    @test abs(expected_nu - nu)*180/pi < 1e-2

    # going back to mean anomaly
    M_ = tra_to_man(nu, e)
    @test abs(mod(M - M_, 2*pi))*180/pi < 1e-8

    # --- 
    # circular 
    M = 1.0pi

    nu = man_to_tra(M, 1e-13)
    @test M == nu

    nu = man_to_tra(M, 1e-6)
    @test M == nu


end