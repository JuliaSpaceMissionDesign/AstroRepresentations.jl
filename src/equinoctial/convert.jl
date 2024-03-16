
function convert_state(::Type{Equi}, c::Cart{N}, μ::Number, args...) where N 
    return Equi{N}(convert6_cart_to_equi(c, μ))
end

function convert_state(::Type{Cart}, c::Equi{N}, μ::Number, args...) where N 
    return Cart{N}(convert6_equi_to_cart(c, μ)) 
end

function convert_state(::Type{Equi}, c::Coe{N}, args...) where N 
    return Equi{N}(convert6_coe_to_equi(c))
end

function convert_state(::Type{Coe}, c::Equi{N}, args...) where N 
    return Coe{N}(convert6_equi_to_coe(c))
end
       
"""
    convert6_equi_to_cart(equi::AbstractVector{<:Number}, μ::Number)

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
function convert6_equi_to_cart(equi::AbstractVector{<:Number}, μ::Number)

    @fastmath begin
        @inbounds p, f, g, h, k, L = @views(equi[1:6])

        # Compute eccentricity and true anomaly to check state validity in case of hyperbolic trajectory
        ecc = sqrt( f*f + g*g )
        if ecc > 1.0
            ν = mod2pi(L - atan(g, f))
            if abs(ν) > acos(-1/ecc)
                throw(
                    ErrorException("Invalid true longitude for hyperbolic trajectory!")
                )
            end
        end
    end

    h2 = h*h 
    k2 = k*k
    ss = 1.0 + h2 + k2
    sLon, cLon = sincos(L)
    tmp1 = 2.0 * h * k
    tmp2 = h2 - k2
    tmp3 = cLon + f
    tmp4 = sLon + g

    ipx = ((tmp2 + 1.0) * cLon + tmp1 * sLon)/ss
    ipy = (-(tmp2 - 1.0) * sLon + tmp1 * cLon)/ss
    ipz = 2.0*(h * sLon - k * cLon)/ss
    ivx = (tmp1 * tmp3 - (tmp2 + 1.0) * tmp4)/ss
    ivy = (-tmp1 * tmp4 - (tmp2 - 1.0) * tmp3)/ss
    ivz = 2.0*(h*tmp3 + k*tmp4)/ss

    w = 1.0 + f * cLon + g * sLon
    r = p/w
    v = sqrt(μ/p)

    return SVector{6}(r*ipx, r*ipy, r*ipz, v*ivx, v*ivy, v*ivz)

end

"""
    convert6_cart_to_equi(sv::AbstractVector{<:Number}, μ::Number)

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
@fastmath function convert6_cart_to_equi(rv::AbstractVector{<:Number}, μ::Number)

    @inbounds R = SVector{3}(rv[1], rv[2], rv[3])
    @inbounds V = SVector{3}(rv[4], rv[5], rv[6])

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

"""
    convert6_coe_to_equi(coe::AbstractVector{<:Number})

Convert classical orbital elements state vector to Equinoctial keplerian elements.

### Inputs 
- `coe` -- Keplerian elements -- `L, rad`

### Output 
Equinoctial representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
"""
@fastmath function convert6_coe_to_equi(coe::AbstractVector{<:Number})
    @inbounds sma, ecc, inc, ran, aop, ta = @views(coe[1:6]) 

    # Support variables
    cy, cx = sincos(ran + aop)
    sr, cr = sincos(ran)
    ti2 = tan(0.5 * inc)

    p = sma * (1 - ecc*ecc)
    f = ecc * cx 
    g = ecc * cy 
    h = ti2 * cr 
    k = ti2 * sr
    L = mod2pi(ran + aop + ta)
    
    return SVector{6}(p, f, g, h, k, L)
end

"""
    convert6_equi_to_coe(coe::AbstractVector{<:Number})

Convert Equinoctial keplerian elements to classical orbital elements state vector .

### Inputs 
- `equi` -- Equinoctial elements -- `L, rad`

### Output 
Keplerian representation of the state as a `SVector{6}`.

### References 
- Vallado, David A. - *Fundamentals of astrodynamics and applications*, 2013.
"""
function convert6_equi_to_coe(equi::AbstractVector{<:Number})
    @inbounds p, f, g, h, k, L = @views(equi[1:6])
    f2 = f*f 
    g2 = g*g

    a = p/(1 - f2 - g2)
    e = sqrt(f2 + g2)
    i = 2*atan( sqrt(h*h + k*k) )

    if i ≤ 1e-9
        Ω = 0
    else 
        Ω = atan(k, h)
    end

    if e ≤ 1e-12 
        ω = 0 
        ϕ = 0
    else
        ϕ = atan(g, f)
        ω = mod2pi(ϕ - Ω)
    end
    ν = mod2pi(L - ϕ)
    return SVector{6}(a, e, i, Ω, ω, ν)
end
