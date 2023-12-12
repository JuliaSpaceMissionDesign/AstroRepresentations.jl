export tan_to_eca, tan_to_hya, tan_to_man, tan_to_par, 
       eca_to_tan, hya_to_tan, man_to_tan, par_to_tan

"""
    tan_to_eca(tan::Number, ecc::Number)

Convert true anomaly to eccentric anomaly.
"""
function tan_to_eca(ta::Number, ecc::Number)
    return 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(ta / 2))
end

"""
    tan_to_hya(tan::Number, ecc::Number)

Convert true anomaly to hyperbolic anomaly.
"""
function tan_to_hya(ta::Number, ecc::Number)
    return 2 * atanh(sqrt((ecc - 1) / (ecc + 1)) * tan(ta / 2))
end

"""
    tan_to_par(ta::Number)

Convert true anomaly to parabolic anomaly.
"""
function tan_to_par(ta::Number, args...)
    return tan(ta / 2)
end

"""
    eca_to_tan(eca::Number, ecc::Number)

Convert eccentric anomaly to true anomaly.
"""
function eca_to_tan(eca::Number, ecc::Number)
    return mod(2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(eca / 2)), 2 * π)
end

"""
    hya_to_tan(hca::Number, ecc::Number)

Convert hyperbolic anomaly to true anomaly.
"""
function hya_to_tan(hca::Number, ecc::Number)
    return 2 * atan(sqrt((ecc + 1) / (ecc - 1)) * tanh(hca / 2))
end

"""
    par_to_tan(pca::Number)

Convert parabolic anomaly to true anomaly.
"""
function par_to_tan(pca::Number, args...)
    return 2 * atan(pca)
end

"""
    man_to_tan(man::Number, ecc::Number)

Convert mean anomaly to true anomaly.
"""
function man_to_tan(ma::Number, ecc::Number)
    if ecc < 0
        throw(
            ArgumentError("eccentricity must be greater than 0")
        )
    end
    err = 1
    TOL = 1e-12

    if isapprox(ecc, 0.0; atol=1e-8)  # circular
        ta = ma

    elseif isapprox(abs(ecc - 1), 0.0; atol=1e-8) # parabolic
        s = (π / 2 - atan(3 * ma / 2)) / 2
        w = atan(tan(s)^(1 / 3))
        B = 2 / tan(2 * w)
        ta = par_to_tan(B)

    elseif 0.0 < ecc < 1.0  # elliptic
        if (-π < ma < 0) || ma > π
            E = ma - ecc
        else
            E = ma + ecc
        end

        it = 1
        while err >= TOL && it <= 10
            f = E - ecc * sin(E) - ma
            g = 1 - ecc * cos(E)
            s = f / g
            E = E - s
            err = abs(s)
        end
        ta = eca_to_tan(E, ecc)

    else  # hyperbolic 
        if ecc < 1.6
            if (-π < ma < 0) || ma > π
                H = ma - ecc
            else
                H = ma + ecc
            end
        else
            if ecc < 3.6 && abs(ma) > π
                H = ma - ecc * sign(ma)
            else
                H = ma / (ecc - 1)
            end
        end

        it = 1
        while err >= TOL && it <= 10
            f = ecc * sinh(H) - H - ma
            g = ecc * cosh(H) - 1
            s = f / g
            H = H - s
            err = abs(s)
        end
        ta = hya_to_tan(H, ecc)
    end
    return ta
end

"""
    tan_to_man(tan::Number, ecc::Number)

Convert true anomaly into mean anomaly
"""
function tan_to_man(ta::Number, ecc::Number)
    if isapprox(ecc, 0; atol=1e-8)  # circular 
        ma = ta
    elseif isapprox(abs(ecc - 1), 0.0; atol=1e-8)  # parabolic
        ma = 1 / 2 * tan(ta / 2) + 1 / 6 * tan(ta / 2)^3
    elseif 0.0 < ecc < 1.0  # elliptic
        E = tan_to_eca(ta, ecc)
        ma = E - ecc * sin(E)
    else  # hyperbolic 
        H = tan_to_hya(ta, ecc)
        ma = ecc * sinh(H) - H
    end
    return ma
end