using AstroRepresentations
using Test
using ForwardDiff
using StaticArrays
using LinearAlgebra

include("src/testutil.jl")

@testset "Representations.jl" verbose=true begin
    include("src/angles.jl")
    include("src/keplerian.jl")
    include("src/equinoctial.jl")
end;