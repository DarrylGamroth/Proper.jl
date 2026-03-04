"""Cubic convolution placeholder using bilinear sample fallback."""
function libcconv(a::AbstractMatrix, y::Real, x::Real)
    return bilinear_sample(a, y, x)
end
