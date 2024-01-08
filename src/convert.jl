export convert_state 

function convert_state(::Type{R}, c::R, args...) where {R <: AbstractStateRepr}
    return c
end

function convert_state(::Type{R1}, c::R2, 
    args...) where {R1<:AbstractStateRepr, R2<:AbstractStateRepr}
    throw(
        ErrorException("State conversion not implemented!")
    )
end 
