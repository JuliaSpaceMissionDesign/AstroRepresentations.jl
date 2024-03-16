
function convert6_cart_to_bplane(x, μ)
    cart = Cart(x)
    pos = position(cart)
    vel = velocity(cart)

    r = norm(pos)
    v = norm(vel)
    pos_dir = pos / r

    # V-infinity magnitude and semi-major axis
    v2 = v^2 
    c3 = v2 - 2μ/r
    @assert c3 > 0 "Cartesian state does not represent an hyperbolic trajectory"
    v∞ = sqrt(c3)
    sma = -μ/v∞

    # Momentum
    h = cross(pos, vel)
    h_dir = unitvec(h)

    # Eccentricity vector 
    ecv = (v2 * pos - dot(pos, vel) * vel)/μ - pos_dir
    ecc = norm(ecv)
    ecv_dir = ecv/ecc

    # Minor axis direction, towards +90 deg true anomaly
    ma_dir = cross(h_dir, ecv_dir)

    # Incoming vinf direction
    tan_inf = asymptote_tra(ecc)
    in_vinf_dir = -cos(tan_inf) * ecv_dir + sin(tan_inf) * ma_dir
    pol_in = convert3_cart_to_sphe(in_vinf_dir)

    # Bplane definition 
    z_dir = SVector{3}(0., 0., 1.)
    t_dir = unitvec( cross(in_vinf_dir, z_dir) )
    r_dir = cross(in_vinf_dir, t_dir)

    # Bplane impact vector 
    b = -sma * sqrt( ecc^2 - 1 )
    bvec = b * cross(in_vinf_dir, h_dir)
    bt = dot(t_dir, bvec)
    br = dot(r_dir, bvec)

    # True anomaly 
    tra = atan( dot(pos, ma_dir), dot(pos, ecv_dir) )
    return SVector{6}(v∞, pol_in[2], pol_in[3], bt, br, tra)

end