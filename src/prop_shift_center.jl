"""Circularly shift array by floor(n/2) along each axis."""
function prop_shift_center(a::AbstractMatrix)
    sy = size(a, 1) ÷ 2
    sx = size(a, 2) ÷ 2
    return circshift(a, (sy, sx))
end

function prop_shift_center!(out::AbstractMatrix{T}, a::AbstractMatrix) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    return copyto!(out, prop_shift_center(a))
end

function prop_shift_center!(out::StridedMatrix{T}, a::StridedMatrix) where {T}
    size(out) == size(a) || throw(ArgumentError("output size must match input size"))
    return half_shift_copy!(out, a)
end
