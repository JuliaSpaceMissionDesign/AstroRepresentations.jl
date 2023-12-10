function convert6_cart_to_topo_radec(rv::AbstractVector{<:Number}, 
    rvₛ::AbstractVector{<:Number}, args...)

    @inbounds ρ =  SVector{3}(rv[1]-rvₛ[1], rv[2]-rvₛ[2], rv[3]-rvₛ[3])
    @inbounds dρ =  SVector{3}(rv[4]-rvₛ[4], rv[5]-rvₛ[5], rv[6]-rvₛ[6])

    ρxy2 = ρ[1]*ρ[1] + ρ[2]*ρ[2]
    ρxy = sqrt(ρxy2)
    ρn = sqrt(ρxy2 + ρ[3]*ρ[3])

    dρxy2 = dρ[1]*dρ[1] + dρ[2]*dρ[2]

    δₜ = asin(ρ[3]/ρn)
    if ρxy ≉ 0.0 
        αₜ = atan(ρ[2]/ρxy, ρ[1]/ρxy)
    else
        dρxy = sqrt(dρxy2)
        αₜ = atan(dρ[2]/dρxy, dρ[1]/dρxy)
    end

    δρ = (ρ[1]*dρ[1] + ρ[2]*dρ[2] + ρ[3]*dρ[3])/ρn 
    δαₜ = - (dρ[1]*ρ[2] - dρ[2]*ρ[1])/dρxy2
    δδₜ = (dρ[3] - δρ*sin(δₜ))/ρxy

    return SVector{6}(ρn, αₜ, δₜ, δρ, δαₜ, δδₜ)
end

function convert3_cart_to_topo_radec(rv::AbstractVector{<:Number}, 
    rvₛ::AbstractVector{<:Number}, args...)

    @inbounds ρ =  SVector{3}(rv[1]-rvₛ[1], rv[2]-rvₛ[2], rv[3]-rvₛ[3])
    ρxy2 = ρ[1]*ρ[1] + ρ[2]*ρ[2]
    ρxy = sqrt(ρxy2)
    ρn = sqrt(ρxy2 + ρ[3]*ρ[3])
    δₜ = asin(ρ[3]/ρn)
    αₜ = atan(ρ[2]/ρxy, ρ[1]/ρxy)

    return SVector{3}(ρn, αₜ, δₜ)
end