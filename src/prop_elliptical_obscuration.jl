"""Multiply the current wavefront by an elliptical obscuration."""
function prop_elliptical_obscuration(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    opts = EllipseOptions(real(eltype(wf.field)), kwargs)
    return _apply_shifted_ellipse!(wf, xradius, yradius, xc, yc, opts, true)
end
