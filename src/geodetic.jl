

@fastmath function convert3_carth_to_geod(pos::AbstractVector, 
    R::Number, f::Number, args...; ϵ::Number = 1e-12)

    @inbounds x, y, z = pos[1], pos[2], pos[3]
    sz = sign(z)

    rδs = x*x + y*y

    # Get eccentricity from flattening  
    e² = f * (2 - f)
    sϕ² = z^2 / (rδs + z^2)

    zₙ, cₙ = z, 1.0
    err = ϵ + 1
    while err > ϵ
        c = cₙ

        z² = zₙ^2
        sϕ² = z² / (rδs + z²)

        cₙ = R * e² * sqrt(sϕ² / (1 - e² * sϕ²))
        zₙ = z + cₙ * sz

        err = abs(cₙ - c)
    end

    λ = atan(y, x)
    ϕ = atan(zₙ, sqrt(rδs))
    h = sqrt(rδs + zₙ^2) - R / sqrt(1 - e² * sϕ²)

    return SVector{3}(h, λ, ϕ)
end

@fastmath function convert3_geod_to_carth(gd::AbstractVector{<:Number}, 
    R::Number, f::Number, args...)
    @inbounds h, λ, ϕ = @views(gd[1:3])
    
    e² = (2 - f) * f

    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)

    d = R / sqrt(1 - e² * sϕ^2)
    c = (d + h) * cϕ
    s = (1 - e²) * d

    return SVector{3}(c * cλ, c * sλ, (s + h) * sϕ)
end