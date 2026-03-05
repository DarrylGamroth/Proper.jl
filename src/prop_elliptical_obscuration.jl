function prop_elliptical_obscuration(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    ny, nx = size(wf.field)
    m = ensure_mask_buffer!(wf.workspace.mask, ny, nx)
    prop_ellipse!(m, wf, xradius, yradius, xc, yc; kwargs...)
    _apply_shifted_mask!(wf.field, m; invert=true)
    return wf
end
