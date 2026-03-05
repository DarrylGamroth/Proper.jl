using Serialization

"""Load serialized state and copy all wavefront fields into `wf`."""
function prop_state(wf::WaveFront, path::AbstractString)
    loaded = open(path, "r") do io
        deserialize(io)
    end
    loaded isa WaveFront || throw(ArgumentError("State file does not contain WaveFront"))
    size(wf.field) == size(loaded.field) || throw(ArgumentError("State wavefront size mismatch"))

    copyto!(wf.field, loaded.field)
    wf.wavelength_m = loaded.wavelength_m
    wf.sampling_m = loaded.sampling_m
    wf.z_m = loaded.z_m
    wf.beam_diameter_m = loaded.beam_diameter_m
    wf.z_w0_m = loaded.z_w0_m
    wf.w0_m = loaded.w0_m
    wf.z_rayleigh_m = loaded.z_rayleigh_m
    wf.current_fratio = loaded.current_fratio
    wf.reference_surface = loaded.reference_surface
    wf.beam_type_old = loaded.beam_type_old
    wf.propagator_type = loaded.propagator_type
    wf.rayleigh_factor = loaded.rayleigh_factor
    return wf
end
