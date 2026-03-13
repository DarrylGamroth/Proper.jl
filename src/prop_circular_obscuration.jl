function prop_circular_obscuration(wf::WaveFront, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    opts = EllipseOptions(real(eltype(wf.field)), kwargs)
    return _apply_shifted_ellipse!(wf, radius, radius, xc, yc, opts, true)
end
