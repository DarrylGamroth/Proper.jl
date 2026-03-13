struct WaveFrontSnapshot{A,T,RS,BT,PT}
    field::A
    wavelength_m::T
    sampling_m::T
    z_m::T
    beam_diameter_m::T
    z_w0_m::T
    w0_m::T
    z_rayleigh_m::T
    current_fratio::T
    reference_surface::RS
    beam_type_old::BT
    propagator_type::PT
    rayleigh_factor::T
end

function capture_wavefront_state(wf::WaveFront)
    field = similar(wf.field)
    copyto!(field, wf.field)
    return WaveFrontSnapshot(
        field,
        wf.wavelength_m,
        wf.sampling_m,
        wf.z_m,
        wf.beam_diameter_m,
        wf.z_w0_m,
        wf.w0_m,
        wf.z_rayleigh_m,
        wf.current_fratio,
        wf.reference_surface,
        wf.beam_type_old,
        wf.propagator_type,
        wf.rayleigh_factor,
    )
end

function restore_wavefront_state!(wf::WaveFront, snap::WaveFrontSnapshot)
    copyto!(wf.field, snap.field)
    wf.wavelength_m = snap.wavelength_m
    wf.sampling_m = snap.sampling_m
    wf.z_m = snap.z_m
    wf.beam_diameter_m = snap.beam_diameter_m
    wf.z_w0_m = snap.z_w0_m
    wf.w0_m = snap.w0_m
    wf.z_rayleigh_m = snap.z_rayleigh_m
    wf.current_fratio = snap.current_fratio
    wf.reference_surface = snap.reference_surface
    wf.beam_type_old = snap.beam_type_old
    wf.propagator_type = snap.propagator_type
    wf.rayleigh_factor = snap.rayleigh_factor
    return wf
end
