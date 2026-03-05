"""Compatibility helper for upstream `libszoom` behavior (damped-sinc zoom)."""
function libszoom(a::AbstractMatrix, mag::Real, nout::Integer=0)
    return prop_szoom(a, mag, nout)
end
