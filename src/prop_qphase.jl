"""Apply quadratic phase for curvature c (meters)."""
function prop_qphase(wf::WaveFront, c::Real)
    c == 0 && return wf
    ny, nx = size(wf.field)
    rsqr = fft_order_rsqr_map(ny, nx, wf.sampling_m)
    wf.field .*= cis.((pi / (wf.wavelength_m * float(c))) .* rsqr)
    return wf
end
