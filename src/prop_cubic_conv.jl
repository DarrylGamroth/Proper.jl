"""Cubic convolution placeholder entrypoint."""
function prop_cubic_conv(a::AbstractMatrix, y::Real, x::Real)
    return libcconv(a, y, x)
end
