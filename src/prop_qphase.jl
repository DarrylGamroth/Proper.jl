"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real)
    c == 0 && return wf
    r = prop_radius(wf)
    phase = @. pi / (wf.wavelength_m * c) * (r^2)
    wf.field .*= cis.(phase)
    return wf
end
