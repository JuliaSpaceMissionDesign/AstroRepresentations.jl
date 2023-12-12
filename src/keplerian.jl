export convert6_cart_to_coe, convert6_coe_to_cart, 
       ∂convert6_cart_to_coe

"""
    convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, [args]...)

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
@fastmath function convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, args...)

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
    ∂convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, [args]...)

Convert cartesian state representation into classical orbital elements. Compute also the full
jacobian of the elements wrt the cartesian state.

### Inputs 
- `sv` -- state vector -- `L, T`
- `μ` -- center gravitational parameter  -- `L³/T²`

### Output
Classical Orbital Elements representation of the state as a `SVector{6}` and its jacobian as 
a `SMatrix{6, 6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*. Vol. 12. 
  Springer Science & Business Media, 2001.
- Pasquale, A. - *Multiple Shooting Optimiser (MSO)*. Technical Note 0001, 2022.
"""
@fastmath function ∂convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number, args...)

    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    K = SVector{3}(0., 0., 1.)
    
    r = norm(R)
    v² = dot(V, V)
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)
    N = cross(K, H) # node axis
    p = h²/μ
    I3 = SMatrix{3,3}(
        1., 0., 0., 
        0., 1., 0., 
        0., 0., 1.
    )

    RoV = R * V'
    RdotV = dot(R, V) 
    RdotN = dot(R, N)

    # --------------------------------------------------------------------------------------
    # ECCENTRICITY
    E = ((v²-μ/r).*R - RdotV.*V) ./ μ # eccentricity vector 
    ecc = norm(E)  
    En = E/ecc
 
    # partials 
    ∂E∂R = (v²-μ/r)*I3 + R*R'/r^3 - V*V'/μ
    ∂E∂V = (2*RoV - RoV')/μ - RdotV/μ*I3
    ∂e∂R = ∂E∂R * En
    ∂e∂V = ∂E∂V' * En

    # --------------------------------------------------------------------------------------
    # SEMIMAJOR AXIS 
    sma = p/(1-ecc^2) # semimajor axis
    Energy = .5*v² - μ/r
    E² = Energy^2

    # partials
    ∂a∂R = μ^2/(2E²*r^3)*R 
    ∂a∂V = μ/(2E²)*V

    # --------------------------------------------------------------------------------------
    # INCLINATION
    inc = acos(H[3]/h) # inclination

    SxR = skew(R) 
    SxV = skew(V)
    SxRt = SxR'
    SxVt = SxV'

    circular = ecc < 1e-12 
    @inbounds equatorial = H[1]*H[1] + H[2]*H[2] < 1e-12  

    if !circular && !equatorial 

        @inbounds h12 = sqrt(H[1]*H[1] + H[2]*H[2])
        ∂i∂H = SVector{3}(
            H[1]*H[3]/(h² * h12), H[2]*H[3]/(h² * h12), -h12/h²
        )

        ∂i∂R = SxV * ∂i∂H
        ∂i∂V = SxRt * ∂i∂H  

        # ----------------------------------------------------------------------------------
        # RIGHT-ASCENTION OF THE ASCENDING NODE 
        ran = mod2pi(atan(N[2], N[1]))

        # partials 
        N12 = h12*h12
        ∂Ω∂h = SVector{3}(
            -H[2]*sign(H[1])/N12, abs(H[1])/N12, 0.
        )

        ∂Ω∂R = SxV * ∂Ω∂h
        ∂Ω∂V = SxRt * ∂Ω∂h

        # ----------------------------------------------------------------------------------
        # TRUE ANOMALY 
        if sma > 0 
            # elliptical orbit 
            e_se = RdotV / sqrt(μ*sma)
            e_ce = r * v²/μ - 1
            ta = eca_to_tan(atan(e_se, e_ce), ecc)
        else 
            # hyperbolic orbit
            e_sh = RdotV / sqrt(-μ*sma)
            e_ch = r * v²/μ - 1
            ta = hya_to_tan(log((e_ch+e_sh)/(e_ch-e_sh))/2, ecc)
        end

        # partials
        Rn = R/r 
        ∂nν∂R = (2*v²/μ*r - 1) * Rn - 2/μ * RdotV * V 
        ∂nν∂V = 2r^2/μ*V - 2/μ * RdotV * R
        ∂dν∂R = ecc*Rn + ∂e∂R*r
        ∂dν∂V = ∂e∂V*r

        fν = dot(En, Rn)
        ∂ν∂f = -1/sqrt(1-fν^2)
        ∂fν∂R = (∂nν∂R - fν*∂dν∂R)/(ecc*r)
        ∂fν∂V = (∂nν∂V - fν*∂dν∂V)/(ecc*r)
        ∂ν∂R = ∂ν∂f * ∂fν∂R
        ∂ν∂V = ∂ν∂f * ∂fν∂V

        # ----------------------------------------------------------------------------------
        # ARGUMENT OF PERIGEE
        px = RdotN
        py = dot(R, cross(H, N))/h
        aop = mod2pi(atan(py, px) - ta)

        # partials
        n = norm(N)
        Nn = N/n
        fω = dot(Nn, En)
        ∂ω∂f = 1/sqrt(1-fω^2)

        SxH = skew(H) 
        SxHt = SxH'
        SxE = skew(E)
        DhN = SMatrix{3, 3}(
            0.,     1.,     0., 
            -1.,    0.,     0., 
            0.,     0.,     0.
        )
        ∂N∂R = DhN * SxVt 
        ∂N∂V = DhN * SxR
        ∂nω∂R = - (SxHt * ∂E∂R + SxE * SxVt)' * K 
        ∂nω∂V = - (SxHt * ∂E∂V + SxE * SxR)' * K
        ∂dω∂R = ∂N∂R' * Nn * ecc + n * ∂e∂R
        ∂dω∂V = ∂N∂V' * Nn * ecc + n * ∂e∂V

        ∂ω∂R = ∂ω∂f * (∂nω∂R - fω*∂dω∂R)/(ecc*n)
        ∂ω∂V = ∂ω∂f * (∂nω∂V - fω*∂dω∂V)/(ecc*n)

    else
        throw(
            ErrorException("orbit is either equatorial or circular: cart_to_coe jacobian ill defined")
        )
    end
    
    kep = SVector{6}(sma, ecc, inc, ran, aop, ta)
    ∂kep∂s = SMatrix{6, 6}(
        ∂a∂R[1],    ∂e∂R[1],    ∂i∂R[1],    ∂Ω∂R[1],     ∂ω∂R[1],     ∂ν∂R[1],
        ∂a∂R[2],    ∂e∂R[2],    ∂i∂R[2],    ∂Ω∂R[2],     ∂ω∂R[2],     ∂ν∂R[2],
        ∂a∂R[3],    ∂e∂R[3],    ∂i∂R[3],    ∂Ω∂R[3],     ∂ω∂R[3],     ∂ν∂R[3],
        ∂a∂V[1],    ∂e∂V[1],    ∂i∂V[1],    ∂Ω∂V[1],     ∂ω∂V[1],     ∂ν∂V[1],
        ∂a∂V[2],    ∂e∂V[2],    ∂i∂V[2],    ∂Ω∂V[2],     ∂ω∂V[2],     ∂ν∂V[2],
        ∂a∂V[3],    ∂e∂V[3],    ∂i∂V[3],    ∂Ω∂V[3],     ∂ω∂V[3],     ∂ν∂V[3],
    )
    return kep, ∂kep∂s
end

"""
    convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, [args]...)

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
function convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, args...) 
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

# """
#     ∂convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, [args]...)

# Convert classical orbital elements state vector to cartesian state. Compute also the full
# jacobian of the cartesian states wrt the elements.

# ### Inputs 
# - `sv` -- cartesian state representation -- `L, rad`
# - `μ` -- Center's gravitational parameter  -- `L³/T²`

# ### Output 
# Cartesian representation of the state as a `SVector{6}` and its jacobian as a `SMatrix{6, 6}`.

# ### References 
# - Vallado, David A. - *Fundamentals of astrodynamics and applications*. Vol. 12. 
#   Springer Science & Business Media, 2001.
# - Pasquale, A. - *Multiple Shooting Optimiser (MSO)*. Technical Note 0001, 2022.
# """ # FIXME: not currently working for ran, aop, tan derivatives
# @fastmath function ∂convert6_coe_to_cart(sv::AbstractVector{<:Number}, μ::Number, args...) 
#     @inbounds sma, ecc, inc, ran, aop, ta = @views(sv[1:6])  
#     e² = ecc*ecc

#     ∂p∂e = -2*sma*ecc
#     ∂p∂a = 1-e²
    
#     p = sma*∂p∂a

#     stan, ctan = sincos(ta)
#     r = p/(1 + ecc*ctan)
#     v = sqrt(2μ/r - μ/sma)
#     sran, cran = sincos(ran)
#     sinc, cinc = sincos(inc)
#     saop, caop = sincos(aop)
#     s1 = saop*cinc 
#     s2 = caop*cinc

#     c1x = caop*cran - s1*sran
#     c1y = saop*cran + s2*sran
#     c2x = caop*sran + s1*cran 
#     c2y = s2*cran - saop*sran
#     c3x = saop*sinc 
#     c3y = caop*sinc

#     rx = r*ctan
#     ry = r*stan
#     μop = sqrt(μ/p) 
#     vx = -μop*stan
#     vy = μop*(ecc + ctan)

#     R = SVector{3}(rx*c1x - ry*c1y, rx*c2x + ry*c2y, rx*c3x + ry*c3y) 
#     V = SVector{3}(vx*c1x - vy*c1y, vx*c2x + vy*c2y, vx*c3x + vy*c3y) 

#     # partials (sma)
#     ∂R∂a = R/sma 
#     ∂V∂a = μ/(2*sma*v^2)*(1/sma-2/r)*V

#     # partials (ecc)
#     tmp = p/r
#     ∂r∂e = (∂p∂e*tmp - ctan*p)/tmp^2
#     ∂R∂e = ∂r∂e * R/r 
#     ∂vx∂e = ecc/(1-e²)*vx
#     ∂vy∂e = ecc/(1-e²)*vy + μop

#     # partials (inc)
#     ∂c1x∂i = c3x*sran
#     ∂c1y∂i = -c3y*sran
#     ∂c2x∂i = -c3x*cran
#     ∂c2y∂i = -c3y*cran
#     ∂c3x∂i = s1
#     ∂c3y∂i = s2

#     # partials (ran)
#     ∂c1x∂Ω = -c2x
#     ∂c1y∂Ω = c2y
#     ∂c2x∂Ω = c1x
#     ∂c2y∂Ω = -c1y

#     # partials (aop)
#     ∂c1x∂ω = -c1y
#     ∂c1y∂ω = c1x
#     ∂c2x∂ω = c2y
#     ∂c2y∂ω = -c2x
#     ∂c3x∂ω = c3y
#     ∂c3y∂ω = -c3x

#     # partials (tan)
#     den = tmp*tmp
#     ∂rx∂ν = -p*stan/den
#     ∂ry∂ν = p*(ecc+ctan)/den
#     ∂vx∂ν = -μop*ctan
#     ∂vy∂ν = -μop*stan

#     car = vcat(R, V)
#     ∂car = SMatrix{6, 6}(
#         ∂R∂a[1],  ∂R∂e[1],              ∂c1x∂i*rx-∂c1y∂i*ry,  ∂c1x∂Ω*rx-∂c1y∂Ω*ry,  ∂c1x∂ω*rx-∂c1y∂ω*ry,  ∂rx∂ν*c1x - ∂ry∂ν*c1y,
#         ∂R∂a[2],  ∂R∂e[2],              ∂c2x∂i*rx+∂c2y∂i*ry,  ∂c2x∂Ω*rx+∂c2y∂Ω*ry,  ∂c2x∂ω*rx+∂c2y∂ω*ry,  ∂rx∂ν*c2x + ∂ry∂ν*c2y,
#         ∂R∂a[3],  ∂R∂e[3],              ∂c3x∂i*rx+∂c3y∂i*ry,  0.,                   ∂c3x∂ω*rx+∂c3y∂ω*ry,  ∂rx∂ν*c3x + ∂ry∂ν*c3y,
#         ∂V∂a[1],  ∂vx∂e*c1x-∂vy∂e*c1y,  ∂c1x∂i*vx-∂c1y∂i*vy,  ∂c1x∂Ω*vx-∂c1y∂Ω*vy,  ∂c1x∂ω*vx-∂c1y∂ω*vy,  ∂vx∂ν*c1x - ∂vy∂ν*c1y,
#         ∂V∂a[2],  ∂vx∂e*c2x+∂vy∂e*c2y,  ∂c2x∂i*vx+∂c2y∂i*vy,  ∂c2x∂Ω*vx+∂c2y∂Ω*vy,  ∂c2x∂ω*vx+∂c2y∂ω*vy,  ∂vx∂ν*c2x + ∂vy∂ν*c2y,
#         ∂V∂a[3],  ∂vx∂e*c3x+∂vy∂e*c3y,  ∂c3x∂i*vx+∂c3y∂i*vy,  0.,                   ∂c3x∂ω*vx+∂c3y∂ω*vy,  ∂vx∂ν*c3x + ∂vy∂ν*c3y,
#     )'
#     return car, ∂car
# end


