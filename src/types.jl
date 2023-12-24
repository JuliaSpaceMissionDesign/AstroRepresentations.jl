
export Coe, Cart6, Cart3, Sphe6, Sphe3, Equi

# -----
# Abstract 

abstract type AbstractRepresentationX{D, N} <: FieldVector{D, N} end

abstract type AbstractRepresentation3{N} <: AbstractRepresentationX{3, N} end

abstract type AbstractRepresentation6{N} <: AbstractRepresentationX{6, N} end

# -----
# 6D Representations

struct Coe{N} <: AbstractRepresentation6{N}
    sma::N 
    ecc::N 
    inc::N 
    ran::N
    aop::N 
    tan::N
end

struct Cart6{N} <: AbstractRepresentation6{N}
    pox::N 
    poy::N 
    poz::N 
    vex::N 
    vey::N 
    vez::N
end

struct Sphe6{N} <: AbstractRepresentation6{N}
    r::N 
    ras::N 
    dec::N 
    dr::N 
    dras::N 
    ddec::N
end

struct Equi{N} <: AbstractRepresentation6{N} 
    slr::N 
    ecx::N 
    ecy::N 
    inx::N 
    iny::N 
    lon::N
end
