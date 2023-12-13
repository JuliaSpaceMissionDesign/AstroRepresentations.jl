export convert3_sphe_to_cart, convert3_cart_to_sphe,
       convert6_cart_to_sphe, convert6_sphe_to_cart

"""
    convert6_cart_to_sphe(sv::AbstractVector{<:Number}, [args]...)

Convert cartesian state into spherical (ra-dec) representation.
"""
@fastmath function convert6_cart_to_sphe(sv::AbstractVector{<:Number}, args...) 
    @inbounds px, py, pz, vx, vy, vz = @views(sv[1:6]) 
    rxy2 = px*px + py*py
    rxy = sqrt(rxy2)
    r = sqrt(rxy2 + pz*pz)

    δ = atan(pz/r, rxy/r)
    if rxy ≉ 0.0 
        α = atan(py/rxy, px/rxy)
    else 
        vxy = sqrt(vx*vx + vy*vy)
        α = atan(vy/vxy, vx/vxy)
    end

    dr = (px*vx + py*vy + pz*vz)/r
    dα = -(vx*py - vy*px)/rxy2 
    dδ = (vz - dr*pz/r)/rxy

    return SVector{6}(r, α, δ, dr, dα, dδ)
end

"""
    convert6_sphe_to_cart(sv::AbstractVector{<:Number}, [args]...)

Convert spherical (ra-dec) state into cartesian representation.
"""
@fastmath function convert6_sphe_to_cart(rd::AbstractVector{<:Number}, args...)
    @inbounds r, α, δ, dr, dα, dδ = @views(rd[1:6])

    sα, cα = sincos(α)
    sδ, cδ = sincos(δ)

    cδcα = cδ*cα
    cδsα = cδ*sα

    return SVector{6}(
        r*cδcα,
        r*cδsα,
        r*sδ,
        dr*cδcα - r*sδ*cα*dδ - r*cδsα*dα,
        dr*cδsα - r*sδ*sα*dδ + r*cδcα*dα,
        dr*sδ + r*cδ*dδ
    )
end

"""
    convert3_cart_to_sphe(sv::AbstractVector{<:Number}, [args]...)

Convert cartesian position into spherical (ra-dec).
"""
@fastmath function convert3_cart_to_sphe(sv::AbstractVector{<:Number}, args...) 
    @inbounds px, py, pz = @views(sv[1:3]) 
    rxy2 = px*px + py*py
    rxy = sqrt(rxy2)
    r = sqrt(rxy2 + pz*pz)
    δ = atan(pz/r, rxy/r)
    α = atan(py/rxy, px/rxy)
    return SVector{3}(r, α, δ)
end

"""
    convert6_sphe_to_cart(sv::AbstractVector{<:Number}, [args]...)

Convert spherical (ra-dec) position into cartesian.
"""
@fastmath function convert3_sphe_to_cart(rd::AbstractVector{<:Number}, args...)
    @inbounds r, α, δ = @views(rd[1:3])
    sα, cα = sincos(α)
    sδ, cδ = sincos(δ)

    return SVector{3}(r*cδ*cα, r*cδ*sα, r*sδ)
end