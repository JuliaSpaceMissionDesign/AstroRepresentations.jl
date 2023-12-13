
function randcoe_ell()
    a = rand(0.1:0.1:10.0)
    e = rand(0.01:0.01:0.975)
    i = rand(-π/2:π/18:π/2)
    Ω = rand(-π:π/18:π)
    ω = rand(-π:π/18:π)
    ν = rand(-π:π/18:π)
    return SVector{6}(a, e, i, Ω, ω, ν)
end

function randcart()

    pz = rand(-1:0.01:1)
    α = rand(-2π:π/18:2π)
    d = rand(0.1:0.1:100)
    γ = rand(0.1:0.01:0.98*sqrt(2))

    vC = sqrt(1.0/d) * γ
    px = cos(α) * d 
    py = sin(α) * d 
    vx = -sin(α) * vC
    vy = cos(α) * vC
    vz = rand(-1:0.01:1)
    return SVector{6}(px, py, pz, vx, vy, vz)

end