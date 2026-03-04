"""Normalized sinc: sin(pi*x)/(pi*x)."""
function prop_sinc(x)
    return iszero(x) ? one(float(x)) : sinpi(float(x)) / (pi * float(x))
end
