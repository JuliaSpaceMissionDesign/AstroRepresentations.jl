export convert6_cart_to_equi, convert6_equi_to_cart, 
       convert6_coe_to_equi

"""
    convert6_cart_to_equi(sv::AbstractVector{<:Number}, μ::Number, [args]...)

Convert cartesian state vector to Equinoctial keplerian elements.

### Inputs 
- `sv` -- state vector -- `L, T`
- `μ` -- center gravitational parameter  -- `L³/T²`

### Output 
Equinoctial representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
- Cefola - DOI: 10.2514/6.1972-937
"""
convert6_cart_to_equi

"""
    convert6_equi_to_cart(equi::AbstractVector{<:Number}, μ::Number, [args]...)

Convert equinoctial state elements to cartesian state.

### Inputs 
- `equi` -- Equinoctial elements -- `L, rad`
- `μ` -- center gravitational parameter  -- `L³/T²`

### Output 
Cartesian representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
- Cefola - DOI: 10.2514/6.1972-937
"""
convert6_equi_to_cart

"""
    convert6_coe_to_equi(coe::AbstractVector{<:Number}, [args]...)

Convert classical orbital elements state vector to Equinoctial keplerian elements.

### Inputs 
- `coe` -- Keplerian elements -- `L, rad`

### Output 
Equinoctial representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
"""
convert6_coe_to_equi

"""
    ∂convert6_coe_to_equi(coe::AbstractVector{<:Number}, [args]...)

Convert classical orbital elements state vector to Equinoctial keplerian elements. 
Compute also the full jacobian of the equinoctial elemenents wrt the classical
orbital elements.

### Inputs 
- `coe` -- Keplerian elements -- `L, rad`

### Output 
Equinoctial representation of the state as a `SVector{6}` and its jacobian as 
a `SMatrix{6, 6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
"""
∂convert6_coe_to_equi


function convert6_equi_to_cart(equi::AbstractVector{<:Number}, μ::Number, args...)

    @fastmath begin
        @inbounds slr, ecx, ecy, inx, iny, tlo = @views(equi[1:6])

        # Compute eccentricity and true anomaly to check state validity in case of hyperbolic trajectory
        ecc = sqrt( ecx*ecx + ecy*ecy )
        if ecc > 1.0
            ν = mod2pi(tlo - atan(ecy, ecx))
            if abs(ν) < acos(-1/ecc)
                throw(
                    ErrorException("Invalid true longitude for hyperbolic trajectory!")
                )
            end
        end
    end

    inx2 = inx*inx 
    iny2 = iny*iny
    ss = 1.0 + inx2 + iny2
    sLon, cLon = sincos(tlo)
    tmp1 = 2.0 * inx * iny
    tmp2 = inx2 - iny2
    tmp3 = cLon + ecx
    tmp4 = sLon + ecy

    ipx = ((tmp2 + 1.0) * cLon + tmp1 * sLon)/ss
    ipy = (-(tmp2 - 1.0) * sLon + tmp1 * cLon)/ss
    ipz = 2.0*(inx * sLon - iny * cLon)/ss
    ivx = (tmp1 * tmp3 - (tmp2 + 1.0) * tmp4)/ss
    ivy = (-tmp1 * tmp4 - (tmp2 - 1.0) * tmp3)/ss
    ivz = 2.0*(inx*tmp3 + iny*tmp4)/ss

    w = 1.0 + ecx * cLon + ecy * sLon;
    r = slr/w;
    v = sqrt(μ/slr);

    return SVector{6}(r*ipx, r*ipy, r*ipz, v*ivx, v*ivy, v*ivz)

end

@inbounds @fastmath function convert6_cart_to_equi(rv::AbstractVector{<:Number}, μ::Number, args...)

    R = SVector{3}(rv[1], rv[2], rv[3])
    V = SVector{3}(rv[4], rv[5], rv[6])

    H = cross(R, V)
    rn = norm(R)
    hn = norm(H)
    r̂ = R/rn
    ĥ = H/hn
    
    RdV = dot(R, V)
    v̂ = (rn*V - RdV * r̂)/hn

    p = hn^2/μ
    k = ĥ[1]/(1 + ĥ[3])
    h = -ĥ[2]/(1 + ĥ[3])

    k2 = k*k 
    h2 = h*h
    s2 = 1 + k2 + h2 
    tkh = 2 * k * h
    E = cross(V, H)/μ - r̂
    F = SVector{3}(1 - k2 + h2, tkh, -2*k) / s2
    G = SVector{3}(tkh, 1 + k2 - h2, 2*h) / s2

    f = dot(E, F)
    g = dot(E, G)
    L = atan(r̂[2] - v̂[1], r̂[1] + v̂[2])

    return SVector{6}(p, f, g, h, k, L)

end

@fastmath function convert6_coe_to_equi(coe::AbstractVector{<:Number}, args...)
    @inbounds sma, ecc, inc, ran, aop, ta = @views(coe[1:6]) 

    # Support variables
    cy, cx = sincos(ran + aop)
    sr, cr = sincos(ran)
    ti2 = tan(0.5 * inc)

    slr = sma * (1 - ecc*ecc)
    ecx = ecc * cx 
    ecy = ecc * cy 
    inx = ti2 * cr 
    iny = ti2 * sr
    tlo = mod2pi(ran + aop + ta)
    
    return SVector{6}(slr, ecx, ecy, inx, iny, tlo)
end

@fastmath function ∂convert6_coe_to_equi(coe::AbstractVector{<:Number}, args...)
    @inbounds sma, ecc, inc, ran, aop, ta = @views(coe[1:6]) 

    # Support variables
    cy, cx = sincos(ran + aop)
    sr, cr = sincos(ran)
    ti2 = tan(0.5 * inc)

    slr = sma * (1 - ecc*ecc)
    ecx = ecc * cx 
    ecy = ecc * cy 
    inx = ti2 * cr 
    iny = ti2 * sr
    tlo = mod2pi(ran + aop + ta)

    # Derivatives
    ∂slr_∂sma = 1 - ecc*ecc
    ∂slr_∂ecc = -2*sma*ecc
    ∂ecx_∂ecc = cx
    ∂ecx_∂ran = -ecc*cy
    ∂ecx_∂aop = -ecc*cy
    ∂ecy_∂ecc = cy
    ∂ecy_∂ran = ecc*cx
    ∂ecy_∂aop = ecc*cx
    ∂inx_∂inc = 0.5*(1 + ti2*ti2)*cr
    ∂inx_∂ran = -ti2*sr
    ∂iny_∂inc = 0.5*(1 + ti2*ti2)*sr
    ∂iny_∂ran = ti2*cr
    ∂tlo_∂ran = 1
    ∂tlo_∂aop = 1
    ∂tlo_∂ta = 1
    
    equi = SVector{6}(slr, ecx, ecy, inx, iny, tlo)
    ∂equi = SMatrix{6, 6}(
        ∂slr_∂sma, ∂slr_∂ecc,         0,         0,         0,        0,
                0, ∂ecx_∂ecc,         0, ∂ecx_∂ran, ∂ecx_∂aop,        0,
                0, ∂ecy_∂ecc,         0, ∂ecy_∂ran, ∂ecy_∂aop,        0,
                0,         0, ∂inx_∂inc, ∂inx_∂ran,         0,        0,
                0,         0, ∂iny_∂inc, ∂iny_∂ran,         0,        0,
                0,         0,         0, ∂tlo_∂ran, ∂tlo_∂aop, ∂tlo_∂ta
    )'
    return equi, ∂equi
end
