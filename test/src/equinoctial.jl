@testset "Equinoctial" verbose=true begin

    @testset "conversion" begin
        # Vallado, D. A. (2013), pg. 114
        μ = 3.986004415e5  # km^3/s^2
        rv = SVector{6}(6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341)
        kep = SVector{6}(
            36127.343, 0.832853, 
            deg2rad(87.869126), deg2rad(227.8982603), deg2rad(53.384930), deg2rad(92.3351567))
        
        eq0 = convert6_coe_to_equi(kep)
        eq = convert6_cart_to_equi(rv, μ)

        @test eq[1] ≈ eq0[1] atol=1.0 
        for i in 2:length(eq)
            @test eq[i] ≈ eq0[i] atol=1e-6
        end
        
    end

    @testset "consistency" verbose=true begin
        
        @testset "equatorial/circular orbit" begin
            sv = SVector{6}(cos(π/4), sin(π/4), 0.0, -sin(π/4), cos(π/4), 0.0)
            equi = convert6_cart_to_equi(sv, 1.0)
            eq0 = SVector{6}(1.0, 0.0, 0.0, 0.0, 0.0, π/4)
            for i in eachindex(eq0)
                @test equi[i] ≈ eq0[i] atol=1e-12
            end
        end

        @testset "equatorial orbit" begin
            sv = convert6_coe_to_cart([1.0, 0.05, 0.0, 0.0, 0.0, 0.0], 1.0)
            equi = convert6_cart_to_equi(sv, 1.0)

            @test equi[1] ≈ 1-0.05^2    atol=1e-12
            @test equi[2] ≈ 0.05        atol=1e-12 
            @test equi[3] ≈ 0.0         atol=1e-12
            
            sv = convert6_coe_to_cart([1.0, 0.05, 0.0, π/4, 0.0, 0.0], 1.0)
            equi = convert6_cart_to_equi(sv, 1.0)

            @test sqrt(equi[2]^2 + equi[3]^2) ≈ 0.05   atol=1e-12 
            @test equi[end] ≈ π/4       atol=1e-12
        end

        @testset "circular polar orbit" begin
            sv = SVector{6}(1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
            equi = convert6_cart_to_equi(sv, 1.0)
            eq0 = SVector{6}(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            for i in eachindex(eq0)
                @test equi[i] ≈ eq0[i] atol=1e-12
            end
        end

    end

end
