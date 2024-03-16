
function convert_state(::Type{Vinf}, c::Cart{N}, μ::Number, args...) where N 
    return Vinf{N}(convert6_cart_to_vinf(c, μ, args...))
end

function convert_state(::Type{Cart}, c::Vinf{N}, μ::Number, args...) where N 
    return Cart{N}(convert6_vinf_to_cart(c, μ))
end

function rpe_to_ecc(rp::Number, v∞::Number, μ::Number)
    return 1 + rp * v∞^2/μ
end

function deflection_to_rpe(δ::Number, v∞, μ)
    v∞² = v∞^2
    b = μ / tan(0.5 * δ) / v∞²
    return sqrt( (μ/v∞²)^2 + b^2 ) - μ/v∞²
end

function rpe_to_deflection(rpe, v∞, μ)
    sind = 1/(1 + rpe * v∞^2/μ)
    return 2*abs( asin(sind) )
end

function imp_to_rpe(b, v∞, μ)
    a = -μ/v∞^2
    e = sqrt(1 + (b/a)^2)
    return a*(1 - e)
end

function rpe_to_imp(rpe, v∞, μ)
    a = -μ/v∞^2 
    e = 1 - rpe/a
    return - a * sqrt(e^2 - 1)
end

function rpe_to_vpe(rpe, v∞, μ)
    return sqrt( v∞^2 + 2/μ * rpe)
end

function asymptote_tra(ecc::Number)
    @assert ecc > 1 "eccentricity shall belong to an hyperbolic trajectory"
    return acos( -1.0/ecc )
end

function convert6_cart_to_vinf(cart, μ)

    equi = convert6_cart_to_equi(cart, μ)
    arg = atan(equi[3], equi[2])
    tra = equi[6] - arg
    equi = SVector{6}(equi[1], equi[2], equi[3], equi[4], equi[5], arg)
    cart_peri = convert6_equi_to_cart(equi, μ)

    P = @views cart_peri[1:3]
    V = @views cart_peri[4:6]
    r = norm(P)
    v = norm(V)
    uₚ = P / r
    vₚ = V / v

    v2 = v^2
    vinfSqr = v2 - 2.0 * μ / r

    @assert vinfSqr > 0.0 "Cartesian state does not represent hyperbolic trajectory"

    v∞ = sqrt(vinfSqr)
    ecc = rpe_to_ecc(r, v∞, μ)
    tan_inf = asymptote_tra(ecc)
    ctan = cos(tan_inf)
    stan = sin(tan_inf)

    in_dir = -ctan * uₚ + stan * vₚ
    out_dir = ctan * uₚ + stan * vₚ

    in_pol = convert3_cart_to_sphe(in_dir)
    out_pol = convert3_cart_to_sphe(out_dir) 

    return SVector{6}(in_pol[2], in_pol[3], v∞, out_pol[2], out_pol[3], tra)

end


function convert6_vinf_to_cart(vinf, μ)

    v∞ = vinf[3]
    pol_in = SVector{3}(v∞, vinf[1], vinf[2])
    pol_out = SVector{3}(v∞, vinf[4], vinf[5])

    vinf_in = convert3_sphe_to_cart(pol_in)
    vinf_out = convert3_sphe_to_cart(pol_out)

    # Incoming and outgoing V-infty unit vectors
    in_dir = unitvec(vinf_in)
    out_dir = unitvec(vinf_out)

    # Position and velocity unit vectors at pericenter 
    pos_dir_p = unitvec(in_dir - out_dir)
    vel_dir_p = unitvec(in_dir + out_dir)

    # Position vector at the given true anomaly can be computed from the bases vectors
    stan, ctan = sincos(vinf[6])
    pos_dir = ctan * pos_dir_p + stan * vel_dir_p

    # Periapsis radius and distance magnitude at given true anomaly can be 
    # computed from deflection
    δ = acos( dot(in_dir, out_dir) )
    rpe = deflection_to_rpe(δ, v∞, μ)
    sma = -μ/v∞^2 
    ecc = 1 - rpe/sma
    slr = sma * (1 - ecc^2)
    r = slr/(1 + ecc*ctan )

    # Flight path angle can be used to compute velocity direction at true anomaly 
    # and application of vis-viva yields velocity magnitude
    fpa = -atan( ecc * stan, 1 + ecc * ctan )  + π/2
    stf, ctf = sincos(vinf[6] + fpa)
    vel_dir = ctf * pos_dir_p + stf * vel_dir_p
    v = sqrt( μ * (2/r - 1/sma) )

    pos = r * pos_dir
    vel = v * vel_dir

    return SVector{6}(hcat(pos, vel))

end