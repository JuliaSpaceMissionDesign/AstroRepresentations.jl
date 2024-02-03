@fastmath function ∂convert6_cart_to_coe(sv::AbstractVector{<:Number}, μ::Number)

    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r², v² = dot(R, R), dot(V, V)
    r = sqrt(r²)
    H = cross(R, V) # momentum vector
    h12 = H[1]*H[1] + H[2]*H[2]
    h² = h12 + H[3]*H[3]
    h = sqrt(h²)
    ĥ = H/h
    r̂ = R/r

    E = cross(V, H)/μ - R/r # eccentricity vector 
    p = h²/μ # semilatur rectum
    En = 1/2*v² - μ/r 

    # ---
    # SEMI-MAJOR AXIS
    a = - μ/(2*En) 

    # partials 
    ∂a∂En = - a/En
    ∂a∂R = ∂a∂En * μ * r̂/r²
    ∂a∂V = ∂a∂En * V
    ∂a = vcat(∂a∂R, ∂a∂V)

    # --- 
    # ECCENTRICITY 
    e = norm(E)

    # partials 
    I₃ = SMatrix{3,3}(
        1., 0., 0., 
        0., 1., 0., 
        0., 0., 1.
    )

    SxR = skew(R)
    SxV = skew(V)
    SxH = skew(H)

    ∂r̂∂R = 1/r*(I₃ - r̂*r̂')
    ∂E∂R = -1/μ * SxV * SxV - ∂r̂∂R
    ∂E∂V = -1/μ * (SxH - SxV*SxR)
    ∂E = hcat(∂E∂R, ∂E∂V)

    # partials
    ∂E∂E = E/e
    ∂Ẽ∂E = hcat(∂E∂E, ∂E∂E)
    ∂ẽ = ∂E' * ∂Ẽ∂E
    ∂e = SVector{6}(@views(∂ẽ[:, 1]))

    # --- 
    # INCLINATION
    i = acos(ĥ[3])

    # partials 
    ∂i∂H = - 1/(h²*sqrt(h12)) * SVector{3}(-H[1]*H[3], -H[2]*H[3], h12)
    ∂i = vcat(SxV * ∂i∂H, -SxR * ∂i∂H)

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        rh12 = R[2]*H[1] - R[1]*H[2]
        λ = atan(h*R[3], rh12)
        Ω = mod2pi( atan(ĥ[1], -ĥ[2]) )

        # partials 
        g = h*R[3]/rh12
        rh12² = rh12*rh12
        ∂λ∂g = 1/(1 + g^2) 
        ∂g∂R = h/rh12² * SVector{3}(R[3]*H[2], -R[3]*H[1], rh12)
        ∂g∂H = R[3]/(h*rh12²) * SVector{3}(
            -(H[2]^2*R[2] + H[1]*H[2]*R[1] + H[3]^2*R[2]),
            H[1]^2*R[1] + H[1]*H[2]*R[2] + H[3]^2*R[1],
            H[3] * rh12
        )
        ∂λ∂R = ∂λ∂g * ∂g∂R + ∂λ∂g * SxV * ∂g∂H
        ∂λ∂V = - ∂λ∂g * SxR * ∂g∂H
        ∂λ = vcat(∂λ∂R, ∂λ∂V)

        ∂Ω∂H = 1/h12 * SVector{3}(H[2], -H[1], 0)
        ∂Ω = vcat(-SxV * ∂Ω∂H, SxR * ∂Ω∂H)
    else 
        f1 = atan(R[2], R[1])
        f2 = sign(H[3])
        λ = f1 * f2 # true longitude
        Ω = 0.0

        # partials 
        r12 = R[1]*R[1] + R[2]*R[2]
        ∂f1∂R = f2/r12 * SVector{3}(-R[2], R[1], 0)
        ∂λ = vcat(∂f1∂R, @SVector(zeros(3)))

        ∂Ω = @SVector(zeros(6))
    end

    # check circular orbits 
    if e > 1e-14
        RdV = dot(R, V)
        c = sqrt(p/μ) * 1/(p-r)
        f = atan(sqrt(p/μ)*RdV, p-r)
        ω = mod2pi(λ - f)

        # partials 
        g = c * RdV 
        ∂f∂g = 1/(1 + g^2)
        ∂g∂R = c * V 
        ∂g∂V = c * R
        ∂g∂r = c/(p-r) * RdV 
        ∂g∂p = -c*(r+p)/(2*p*(p-r)) * RdV
        ∂p∂R = -2/μ * SxV * H
        ∂p∂V = 2/μ * SxR * H
        ∂f∂R = ∂f∂g * ( - ∂g∂p * ∂p∂R + ∂g∂r * R/r + ∂g∂R )
        ∂f∂V = ∂f∂g * ( - ∂g∂p * ∂p∂V + ∂g∂V )
        ∂f = vcat(∂f∂R, ∂f∂V)

        ∂ω = ∂λ - ∂f
    else 
        f = mod2pi(λ) 
        ω = 0.

        ∂f = ∂λ
        ∂ω = @SVector(zeros(6))
    end

    ∂M = hcat(∂a, ∂e, ∂i, ∂Ω, ∂ω, ∂f)
    return SA[a, e, i, Ω, ω, f], ∂M'

end
