"""Add wavefront phase error in meters, matching upstream PROPER semantics."""
function prop_add_phase(wf::WaveFront, phase_error::Real)
    wf.field .*= cis.((2pi / wf.wavelength_m) * float(phase_error))
    return wf
end

"""Add wavefront phase-error map in meters at current sampling."""
function prop_add_phase(wf::WaveFront, phase_error::AbstractMatrix)
    size(phase_error) == size(wf.field) || throw(ArgumentError("phase size must match wavefront"))
    if same_backend_style(typeof(wf.field), typeof(phase_error)) && wf.field isa StridedMatrix
        ph = ensure_fft_real_scratch!(wf.workspace.fft, size(phase_error, 1), size(phase_error, 2))
        prop_shift_center!(ph, phase_error; inverse=true)
        @inbounds wf.field .*= cis.((2pi / wf.wavelength_m) .* ph)
    else
        ph = backend_adapt(wf.field, prop_shift_center(phase_error; inverse=true))
        @inbounds wf.field .*= cis.((2pi / wf.wavelength_m) .* ph)
    end
    return wf
end
