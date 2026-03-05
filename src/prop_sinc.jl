"""Sinc interpolation kernel: sin(x)/x with sinc(0)=1."""
function prop_sinc(x)
    xf = float(x)
    return iszero(xf) ? one(xf) : sin(xf) / xf
end
