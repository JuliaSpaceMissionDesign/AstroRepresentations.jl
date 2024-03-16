export convert_state 

function convert_state(::Type{R}, x::R, args...) where {R <: AbstractStateRepr}
    return x
end

function convert_state(::Type{R1}, x::R2, 
    args...) where {R1<:AbstractStateRepr, R2<:AbstractStateRepr}
    throw(ErrorException("State conversion not implemented!"))
end 

function convert_state(t::Type{R}, x::StaticArray{D, N}, 
    args...) where {D, N, R<:AbstractStateRepr6}
    @assert length(x) == 6 "State vector shall be of dimension 6"
    return t(x)
end

function convert_state(t::Type{R}, x::AbstractVector{N}, 
    args...) where {N, R<:AbstractStateRepr6}
    @assert length(x) == 6 "State vector shall be of dimension 6"
    return t(x)
end
