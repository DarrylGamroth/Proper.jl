"""Multiply the current wavefront by a rectangular aperture."""
function prop_rectangular_aperture(wf::WaveFront, width::Real, height::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    ny, nx = size(wf.field)
    m = ensure_mask_buffer!(wf.workspace.mask, ny, nx)
    prop_rectangle!(m, wf, width, height, xc, yc; kwargs...)
    _apply_shifted_mask!(wf.field, m)
    return wf
end
