
function ene(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    v² = dot(V, V)
    r = norm(R)

    return 1/2*v² - μ/r
end

function ∂ene(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r², v² = dot(R, R), dot(V, V)
    r = sqrt(r²)
    r³ = r*r²

    En = 1/2*v² - μ/r

    # partials
    ∂En∂R = μ/r³ * R
    ∂En = vcat(∂En∂R, V)

    return En, ∂En
end

function sma(sv, μ::Number, args...)
    En = ene(sv, μ)
    return - μ/(2*En)
end

function ∂sma(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r², v² = dot(R, R), dot(V, V)
    r = sqrt(r²)
    r³ = r*r²

    En = 1/2*v² - μ/r
    a = - μ/(2*En)

    # partials
    ∂a∂En = - a/En
    ∂a∂R = ∂a∂En * μ * R/r³
    ∂a∂V = ∂a∂En * V
    ∂a = vcat(∂a∂R, ∂a∂V)

    return a, ∂a
end

function mov(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    return cross(R, V)
end

function ∂mov(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V)

    # partials
    ∂H∂R = - skew(V)
    ∂H∂V = skew(R)
    ∂H = hcat(∂H∂R, ∂H∂V)

    return H, ∂H
end

function mom(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V)
    return norm(H)
end

function ∂mom(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V)
    h = norm(H)

    # partials
    SxR = -skew(R)
    SxV = skew(V)

    ∂h∂R = SxV * H / h
    ∂h∂V = SxR * H / h
    ∂h = vcat(∂h∂R, ∂h∂V)

    return h, ∂h
end

function ecv(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    H = cross(R, V) # momentum vector
    E = cross(V, H)/μ - R/r # eccentricity vector 

    return E
end

function ∂ecv(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    ir = R/r
    H = cross(R, V) # momentum vector
    E = cross(V, H)/μ - ir # eccentricity vector 

    I3 = SMatrix{3,3}(
        1., 0., 0., 
        0., 1., 0., 
        0., 0., 1.
    )

    # partials
    SxR = skew(R)
    SxV = skew(V)
    SxH = skew(H)

    ∂ir∂R = 1/r*(I3 - ir*ir')
    ∂E∂R = -1/μ * SxV * SxV - ∂ir∂R
    ∂E∂V = -1/μ * SxH + 1/μ * SxV * SxR
    ∂E = hcat(∂E∂R, ∂E∂V)

    return E, ∂E
end

function ecc(sv, μ::Number, args...)
    return norm(ecv(sv, μ))
end

function ∂ecc(sv, μ::Number, args...)
    E, ∂E = ∂ecv(sv, μ)
    e = norm(E)

    # partials
    ∂E∂E = E/e
    ∂Ẽ∂E = hcat(∂E∂E, ∂E∂E)
    ∂e = ∂E' * ∂Ẽ∂E
    
    return e, SVector{6}(@views(∂e[:, 1]))
end

function inc(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)
    ĥ = H/h
    return acos(ĥ[3]) # inclination 
end

function ∂inc(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    H = cross(R, V) # momentum vector

    h12 = H[1]*H[1] + H[2]*H[2]
    h² = h12 + H[3]*H[3]
    h = sqrt(h²)
    ĥ = H/h
    i = acos(ĥ[3]) # inclination 

    # partials
    ∂H∂R = skew(V)
    ∂H∂V = -skew(R)
    ∂i∂H = - 1/(h²*sqrt(h12)) * SVector{3}(-H[1]*H[3], -H[2]*H[3], h12)
    ∂i = vcat(∂H∂R * ∂i∂H, ∂H∂V * ∂i∂H)

    return i, ∂i
end

function ran(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V) # momentum vector

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        Ω = mod2pi( atan(H[1], -H[2]) )
    else 
        Ω = 0.0
    end

    return Ω
end

function ∂ran(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    H = cross(R, V) # momentum vector
    h12 = H[1]*H[1] + H[2]*H[2]

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        Ω = mod2pi( atan(H[1], -H[2]) )

        # partials
        ∂H∂R = - skew(V)
        ∂H∂V = skew(R)
        ∂Ω∂H = 1/h12 * SVector{3}(H[2], -H[1], 0)

        ∂Ω = vcat(∂H∂R * ∂Ω∂H, ∂H∂V * ∂Ω∂H)
    else 
        Ω = 0.0
        ∂Ω = @SVector(zeros(6))
    end

    return Ω, ∂Ω

end

function aol(sv, args...)

    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        λ = atan(h*R[3], R[2]*H[1] - R[1]*H[2])
    else 
        λ = atan(R[2], R[1]) * sign(H[3]) # true longitude
    end
    return λ

end

function ∂aol(sv, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r12 = R[1]*R[1] + R[2]*R[2]
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)

    # check inclination is zero
    if H[1]*H[1] + H[2]*H[2] > 1e-12 
        rh12 = R[2]*H[1] - R[1]*H[2]
        λ = atan(h*R[3], rh12)

        # partials
        g = h*R[3]/rh12
        rh12² = rh12*rh12

        ∂H∂R = - skew(V)
        ∂H∂V = skew(R)

        ∂λ∂g = 1/(1 + g^2) 
        ∂g∂R = h/rh12² * SVector{3}(R[3]*H[2], -R[3]*H[1], rh12)
        ∂g∂H = R[3]/(h*rh12²) * SVector{3}(
            -(H[2]^2*R[2] + H[1]*H[2]*R[1] + H[3]^2*R[2]),
            H[1]^2*R[1] + H[1]*H[2]*R[2] + H[3]^2*R[1],
            H[3] * rh12
        )
        ∂λ∂R = ∂λ∂g * ∂g∂R - ∂λ∂g * ∂H∂R * ∂g∂H
        ∂λ∂V = - ∂λ∂g * ∂H∂V * ∂g∂H
        ∂λ = vcat(∂λ∂R, ∂λ∂V)

    else 
        f1 = atan(R[2], R[1])
        f2 = sign(H[3])
        λ = f1 * f2 # true longitude

        # partials
        ∂f1∂R = f2/r12 * SVector{3}(-R[2], R[1], 0)
        ∂λ = vcat(∂f1∂R, @SVector(zeros(3)))

    end

    return λ, ∂λ

end

function tra(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    h = sqrt(h²)

    E = cross(V, H)/μ - R/r # eccentricity vector 
    p = h²/μ # semilatur rectum

    ecc = norm(E)
    
    if ecc > 1e-14
        ta = atan(sqrt(p/μ)*dot(R, V), p-r)
    else 
        if H[1]*H[1] + H[2]*H[2] > 1e-12 
            aol = atan(h*R[3], R[2]*H[1] - R[1]*H[2])
        else 
            aol = atan(R[2], R[1]) * sign(H[3]) 
        end
        ta = mod2pi(aol)
    end
    return ta
end

function ∂tra(sv, μ::Number, args...)

    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    H = cross(R, V) # momentum vector
    h² = dot(H, H)

    E = cross(V, H)/μ - R/r # eccentricity vector 
    p = h²/μ # semilatur rectum

    ecc = norm(E)

    if ecc > 1e-14
        RdV = dot(R, V)
        c = sqrt(p/μ) * 1/(p-r)

        ta = atan(sqrt(p/μ)*RdV, p-r)

        g = c * RdV 
        ∂f∂g = 1/(1 + g^2)
        ∂g∂R = c * V 
        ∂g∂V = c * R
        ∂g∂r = c/(p-r) * RdV 
        ∂g∂p = -c*(r+p)/(2*p*(p-r)) * RdV

        ∂H∂R = - skew(V)
        ∂H∂V = skew(R)

        ∂p∂R = 2/μ * ∂H∂R * H
        ∂p∂V = 2/μ * ∂H∂V * H
        
        ∂f∂R = ∂f∂g * ( - ∂g∂p * ∂p∂R + ∂g∂r * R/r + ∂g∂R )
        ∂f∂V = ∂f∂g * ( - ∂g∂p * ∂p∂V + ∂g∂V )
        ∂ta = vcat(∂f∂R, ∂f∂V)
        
    else 
        aol, ∂ta = ∂aol(sv)
        ta = mod2pi(aol)
    end

    return ta, ∂ta
end

function aop(sv, μ::Number, args...) 
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])

    r = norm(R)
    H = cross(R, V) # momentum vector
    h² = dot(H, H)

    E = cross(V, H)/μ - R/r # eccentricity vector 
    e = norm(E)  # eccentricity

    p = h²/μ # semilatur rectum

    if e > 1e-14
        ta = atan(sqrt(p/μ)*dot(R, V), p-r)
        ω = mod2pi(aol(sv) - ta)
    else 
        ω = 0.
    end

    return ω

end

function ∂aop(sv, μ::Number, args...)
    e = ecc(sv, μ)

    if e > 1e-14 

        λ, ∂λ = ∂aol(sv)
        f, ∂f = ∂tra(sv, μ)
        ω = mod2pi(λ - f)
        ∂ω = ∂λ - ∂f

    else 
        ω = 0 
        ∂ω = @SVector(zeros(6))
    end

    return ω, ∂ω
end


function slr(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    p = h²/μ # semilatur rectum
    return p 
end

function ∂slr(sv, μ::Number, args...)
    @inbounds R = SVector{3}(sv[1], sv[2], sv[3])
    @inbounds V = SVector{3}(sv[4], sv[5], sv[6])
    H = cross(R, V) # momentum vector
    h² = dot(H, H)
    p = h²/μ # semilatur rectum

    #partials 
    ∂H∂R = -skew(V)
    ∂H∂V = skew(R)
    ∂p∂R = -2/μ * ∂H∂R * H
    ∂p∂V = -2/μ * ∂H∂V * H
    ∂p = vcat(∂p∂R, ∂p∂V)

    return p, ∂p
end
