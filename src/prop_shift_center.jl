"""Circularly shift array by floor(n/2) along each axis."""
function prop_shift_center(a::AbstractMatrix)
    sy = size(a, 1) ÷ 2
    sx = size(a, 2) ÷ 2
    return circshift(a, (sy, sx))
end
