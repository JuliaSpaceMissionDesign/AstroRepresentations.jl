export convert6_cart_to_coe, convert6_coe_to_cart

"""
    convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number)

Convert cartesian state representation into classical orbital elements.

### Inputs 
- `sv` -- state vector -- `L, T`
- `μ` -- center gravitational parameter  -- `L³/T²`

### Output
Classical orbital elements representation of the state as a `SVector{6}`. 

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*. Vol. 12. 
  Springer Science & Business Media, 2001.
"""
@fastmath function convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number)

    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    v² = dot(V, V)
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)
    ĥ = H/h

    E = cross(V, H)/μ - R/r # eccentricity vector 
    p = h²/μ # semilatur rectum

    sma = (r*μ)/(2*μ - r*v²) # semimajor axis
    ecc = norm(E)  # eccentricity
    inc = acos(ĥ[3]) # inclination 

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        aol = atan(h*R[3], R[2]*H[1] - R[1]*H[2])
        ran = mod2pi( atan(ĥ[1], -ĥ[2]) )
    else 
        aol = atan(R[2], R[1]) * sign(H[3]) # true longitude
        ran = 0.0
    end

    # check circular orbits 
    if ecc > 1e-14
        ta = atan(sqrt(p/μ)*dot(R, V), p-r)
        aop = mod2pi(aol - ta)
    else 
        ta = mod2pi(aol) 
        aop = 0.
    end

    return SVector{6}(sma, ecc, inc, ran, aop, ta)
end

"""
    convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number)

Convert classical orbital elements state vector to cartesian state.

### Inputs 
- `sv` -- Keplerian elements -- `L, rad`
- `μ` -- Center's gravitational parameter  -- `L³/T²`

### Output 
Cartesian representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*. Vol. 12. 
  Springer Science & Business Media, 2001.
"""
function convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number) 
    @inbounds sma, ecc, inc, ran, aop, ta = @views(sv[1:6]) 
    p = sma*(1-ecc*ecc)

    stan, ctan = sincos(ta)
    r = p/(1 + ecc*ctan)
    sran, cran = sincos(ran)
    sinc, cinc = sincos(inc)
    saop, caop = sincos(aop)
    s1 = saop*cinc 
    s2 = caop*cinc

    c1x = caop*cran - s1*sran
    c1y = saop*cran + s2*sran
    c2x = caop*sran + s1*cran 
    c2y = s2*cran - saop*sran
    c3x = saop*sinc 
    c3y = caop*sinc

    rx = r*ctan
    ry = r*stan
    μop = sqrt(μ/p) 
    vx = -μop*stan
    vy = μop*(ecc + ctan)

    return SVector{6}(
        rx*c1x - ry*c1y, rx*c2x + ry*c2y, rx*c3x + ry*c3y,
        vx*c1x - vy*c1y, vx*c2x + vy*c2y, vx*c3x + vy*c3y
    )
end
