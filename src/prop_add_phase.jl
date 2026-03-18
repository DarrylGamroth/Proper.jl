"""Add wavefront phase error in meters, matching upstream PROPER semantics."""
function prop_add_phase(wf::WaveFront, phase_error::Real)
    wf.field .*= cis.((2pi / wf.wavelength_m) * float(phase_error))
    return wf
end

"""Add wavefront phase-error map in meters at current sampling."""
function prop_add_phase(wf::WaveFront, phase_error::AbstractMatrix)
    size(phase_error) == size(wf.field) || throw(ArgumentError("phase size must match wavefront"))
    ph = shift_center_for_wavefront!(wf, phase_error; inverse=true)
    @inbounds wf.field .*= cis.((2pi / wf.wavelength_m) .* ph)
    return wf
end
