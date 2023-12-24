
# Classical orbital elements

function Coe(c::Cart6{N}, μ::Number, args...) where N 
    return Coe{N}(convert6_cart_to_coe(c, μ))
end

function Cart6(c::Coe{N}, μ::Number, args...) where N 
    return Cart6{N}(convert6_coe_to_cart(c, μ))
end

# Equinoctial orbital elements 

function Equi(c::Cart6{N}, μ::Number, args...) where N 
    return Equi{N}(convert6_cart_to_equi(c, μ))
end

function Cart6(e::Equi{N}, μ::Number, args...) where N 
    return Cart6{N}(convert6_equi_to_cart(e, μ))
end

function Equi(c::Coe{N}, args...) where N 
    return Equi{N}(convert6_coe_to_equi(c))
end

# Spherical (ra-dec)

function Sphe6(c::Cart6{N}, args...) where N 
    return Sphe6{N}(convert6_cart_to_sphe(c))
end

function Cart6(c::Sphe6{N}, args...) where N 
    return Cart6{N}(convert6_sphe_to_cart(c))
end