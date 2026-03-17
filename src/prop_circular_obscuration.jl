"""Multiply the current wavefront by a circular obscuration (dark circular mask)."""
function prop_circular_obscuration(wf::WaveFront, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    opts = CircleOptions(real(eltype(wf.field)), kwargs)
    return _apply_shifted_circle!(wf, radius, xc, yc, opts, true)
end
