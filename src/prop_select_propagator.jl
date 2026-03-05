"""Select propagation regime based on Rayleigh-distance boundary."""
function prop_select_propagator(wf::WaveFront, dz::Real)
    dzw = wf.z_w0_m - wf.z_m
    newz = wf.z_m + float(dz)

    beam_type_new = abs(wf.z_w0_m - newz) < wf.rayleigh_factor * wf.z_rayleigh_m ? :INSIDE_ : :OUTSIDE
    wf.propagator_type = Symbol(String(wf.beam_type_old) * "_to_" * String(beam_type_new))
    wf.beam_type_old = beam_type_new
    return dzw
end
