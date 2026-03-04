"""Zoom image with interpolation fallback."""
prop_szoom(a::AbstractMatrix, mag::Real, nout::Integer=0; kwargs...) = prop_magnify(a, mag, nout; kwargs...)
