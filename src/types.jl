
export Coe, Cart, Sphe, Equi

# -----
# Abstract 

abstract type AbstractStateRepr{D, N} <: FieldVector{D, N} end

abstract type AbstractStateRepr6{N} <: AbstractStateRepr{6, N} end

# -----
# 6D Representations

struct Coe{N} <: AbstractStateRepr6{N}
    sma::N 
    ecc::N 
    inc::N 
    ran::N
    aop::N 
    tra::N
end

struct Cart{N} <: AbstractStateRepr6{N}
    pox::N 
    poy::N 
    poz::N 
    vex::N 
    vey::N 
    vez::N
end

struct Sphe{N} <: AbstractStateRepr6{N}
    r::N 
    ras::N 
    dec::N 
    dr::N 
    dras::N 
    ddec::N
end

struct Equi{N} <: AbstractStateRepr6{N} 
    slr::N 
    ecx::N 
    ecy::N 
    inx::N 
    iny::N 
    lon::N
end

function Base.show(io::IO, ::MIME"text/plain", sv::R) where {R<:AbstractStateRepr}
    D = length(sv)
    println(io, "$D-element $(typeof(sv)) with indices SOneTo($D)")
    for el in fieldnames(R)
        val = getfield(sv, el)
        println(io, " :$el => $(val)")
    end
end
