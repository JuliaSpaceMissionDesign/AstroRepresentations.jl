using StaticArrays
struct Point3D <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end
