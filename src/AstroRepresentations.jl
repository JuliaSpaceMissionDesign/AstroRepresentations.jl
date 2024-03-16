module AstroRepresentations

using StaticArrays
using LinearAlgebra
using JSMDUtils.Math: skew, unitvec

include("angles.jl")
include("types.jl")
include("convert.jl")

# Classical orbital elements
include("keplerian/convert.jl")
include("keplerian/jac.jl")
include("keplerian/util.jl")

# Equinoctial orbital elements
include("equinoctial/convert.jl")

# Spherical state representation
include("spherical/convert.jl")

# Vinf state representation
include("vinf.jl")


end
