
export Coe, Cart, Sphe, Equi, Vinf

# -----
# Abstract 

abstract type AbstractStateRepr{D, N <: Number} <: FieldVector{D, N} end

abstract type AbstractStateRepr6{N} <: AbstractStateRepr{6, N} end

# -----
# 6D Representations

"""
    Cart{N} 

Represents a state representation in Cartesian coordinates.
It stores the six Cartesian elements:

- `pox`: Position x-coordinate
- `poy`: Position y-coordinate
- `poz`: Position z-coordinate
- `vex`: Velocity x-coordinate
- `vey`: Velocity y-coordinate
- `vez`: Velocity z-coordinate
"""
struct Cart{N} <: AbstractStateRepr6{N}
    pox::N 
    poy::N 
    poz::N 
    vex::N 
    vey::N 
    vez::N
end

"""
    Coe{N}

Represents a state representation in Classical Orbital Elements (COE).
It stores the six classical orbital elements:

- `sma`: Semi-major axis 
- `ecc`: Eccentricity 
- `inc`: Inclination 
- `ran`: Right ascension of the ascending node 
- `aop`: Argument of periapsis 
- `tra`: True anomaly 
"""
struct Coe{N} <: AbstractStateRepr6{N}
    sma::N 
    ecc::N 
    inc::N 
    ran::N
    aop::N 
    tra::N
end

"""
    Sphe{N}

Represents a state representation in spherical coordinates.
It stores the six spherical elements:

- `rad`: Radius (distance from origin)
- `ras`: Right ascension (RA) or azimuth angle
- `dec`: Declination (DEC) or elevation angle
- `drad`: Rate of change of radius
- `dras`: Rate of change of right ascension or azimuth angle
- `ddec`: Rate of change of declination or elevation angle
"""
struct Sphe{N} <: AbstractStateRepr6{N}
    rad::N 
    ras::N 
    dec::N 
    drad::N 
    dras::N 
    ddec::N
end

"""
    Equi{N}

Represents a state representation in equinoctial elements.
It stores the six equinoctial elements:

- `slr`: Semi-latus rectum 
- `ecx`: Eccentricity vector component in x-direction 
- `ecy`: Eccentricity vector component in y-direction 
- `inx`: Inclination vector component in x-direction 
- `iny`: Inclination vector component in y-direction 
- `lon`: True longitude 
"""
struct Equi{N} <: AbstractStateRepr6{N} 
    slr::N
    ecx::N 
    ecy::N 
    inx::N 
    iny::N 
    lon::N
end

"""
    Vinf{N}

Represents a state representation in V-infinity polar formulation.
It stores the six components:

- `vir`: Incoming V-infinity right-ascesion
- `vid`: Incoming V-infinity declination
- `vinf`: V-infinity magnitue
- `vor`: Outgoing V-infinity right-ascesion
- `vod`: Outgoing V-infinity declination
- `tra`: True anomaly 
"""
struct Vinf{N} <: AbstractStateRepr6{N}
    vir::N 
    vid::N 
    vinf::N 
    vor::N
    vod::N 
    tra::N
end

"""
    BPlane{N}

Represents a state representation in V-infinity b-plane formulation.
It stores the six components:

- `vinf`: V-infinity magnitue
- `vir`: Incoming V-infinity right-ascesion
- `vid`: Incoming V-infinity declination
- `bt`: B-plane T component
- `br`: B-plane R component
- `tra`: True anomaly 
"""
struct BPlane{N} <: AbstractStateRepr6{N}
    vinf::N
    vir::N 
    vid::N 
    bt::N 
    br::N 
    tra::N
end

function Base.show(io::IO, ::MIME"text/plain", sv::R) where {R<:AbstractStateRepr}
    D = length(sv)
    println(io, "$D-element $(typeof(sv)) with indices SOneTo($D)")
    for el in fieldnames(R)
        val = getfield(sv, el)
        println(io, " $el = $(val)")
    end
end

@inline position(x::Cart) = SVector{3}(x[1], x[2], x[3])
@inline velocity(x::Cart) = SVector{3}(x[4], x[5], x[6])
@inline momentum(x::Cart) = cross(position(x), velocity(x))