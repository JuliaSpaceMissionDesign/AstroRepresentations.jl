
@testset "Keplerian" verbose=true begin

    @testset "conversion" begin
        # Vallado, D. A. (2013), pg. 114
        μ = 3.986004415e5  # km^3/s^2
        rv = SVector{6}(6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341)
        kep = SVector{6}(36127.343, 0.832853, 87.869126, 227.8982603, 53.384930, 92.3351567)
        coe = AstroRepresentations.convert6_cart_to_coe(rv, μ)
        rv_ = AstroRepresentations.convert6_coe_to_cart(coe, μ)

        @test all((rv[1:3] - rv_[1:3]) .≤ 1e-8)  # 0.01 mm error
        @test all((rv[3:6] - rv_[3:6]) .≤ 1e-8) # 0.01 mm/s error

        @test isapprox(kep[1], coe[1], atol=1e-2)
        @test isapprox(kep[2], coe[2], atol=1e-6)
        for i in 3:6
            @test isapprox(kep[i], rad2deg(coe[i]), atol=1e-6)
        end
    end
    
    @testset "Consistency" verbose=true begin
        @testset "cart -> coe -> cart" begin
            for μ in (1e3, 1e5, 1e11)
                for pz in (-1., 0., 1.)
                    for α in LinRange(0., 2π, 18)
                        dist = 1.
                        for i in 1:7
                            dist *= 10.
                            for γ in (0.1, 1., 0.99*sqrt(2) )
                                velC = sqrt(μ/dist) * γ
                                px = cos(α) * dist 
                                py = sin(α) * dist 
                                vx = -sin(α) * velC
                                vy = cos(α) * velC
                                
                                cart = SA[px, py, pz, vx, vy, 0.]
                                kep = AstroRepresentations.convert6_cart_to_coe(cart, μ)
                                out = AstroRepresentations.convert6_coe_to_cart(kep, μ)
                                @test all( isapprox.(cart[1:3], out[1:3], atol=1e-8*dist) )
                                @test all( isapprox.(cart[4:end], out[4:end], atol=1e-8) )
                                @test isapprox(norm(cart[1:3]), norm(out[1:3]), atol=1e-10*dist)
                                @test isapprox(norm(cart[4:end]), norm(out[4:end]), atol=1e-8)
                            end
                        end
                    end
                end
            end
        end
    
        @testset "coe -> cart -> coe" begin
            sma = 1.0
            for _ in 1:7
               for ecc in LinRange(0., 0.9, 3)
                    for inc in LinRange(0, π, 6)
                        for ran in LinRange(-π, π, 6)
                            for aop in LinRange(-π, π, 6)
                                for ta in LinRange(-π, π, 6)
                                    coe = SA[sma, ecc, inc, ran, aop, ta]
                                    cart = AstroRepresentations.convert6_coe_to_cart(coe, 1.)
                                    out = AstroRepresentations.convert6_cart_to_coe(cart, 1.)
    
                                    @test isapprox(out[1], coe[1], atol=1e-12)
                                    @test isapprox(out[2], coe[2], atol=1e-9)
                                    @test isapprox(out[3], coe[3], atol=1e-9)
                                    @test (out[4] - mod(coe[4], 2π)) < 1e-12
                                end
                            end
                        end
                    end
               end
            end
        end
    end

    @testset "Jacobians" verbose=true begin 
        @testset "cart->coe" verbose=true begin
            for _ in 1:100
                sv = randcart()
                ref = ForwardDiff.jacobian(
                    x->AstroRepresentations.convert6_cart_to_coe(x, 1.0), sv
                )
                if !any(isnan.(ref))
                    _, jac = AstroRepresentations.∂convert6_cart_to_coe(sv, 1.0)
        
                    if !any(isnan.(jac))
                        @test all(isapprox.(ref, jac; atol=1e-12, rtol=1e-12))
                    else
                        @test_broken all(isapprox.(ref, jac; atol=1e-12, rtol=1e-12))
                    end
                end
            end
        end

    end

    @testset "Utils" verbose=true begin

        funs = [:sma, :ecc, :inc, :ran, :aop, :tra, :slr, :ene, :mom, :aol]
        dfuns = [Symbol("∂$f") for f in funs]
    
        for (f, df) in zip(funs, dfuns)
    
            @testset "$f/$df" verbose=true begin 
                for _ in 1:10
                    @eval begin 
                        # sv = convert_state(Cart, Coe(randcoe_ell()), 1.0)
                        sv = Cart(randcart())
                        valref = AstroRepresentations.$(f)(sv, 1.0)
                        val, ∂val = AstroRepresentations.$(df)(sv, 1.0) 

                        @test isapprox(val, valref; rtol=1e-12)

                        ∂ref = ForwardDiff.gradient(
                            x->AstroRepresentations.$(f)(x, 1.0), 
                            sv
                        )

                        if !any(isnan.(∂ref))
                            @test all(isapprox.(∂ref, ∂val; rtol=1e-8))
                        end
                        
                    end
                end
            end
        end
    end;

end;
