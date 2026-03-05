using Serialization

"""Load serialized state and copy all wavefront fields into `wf`."""
function prop_state(wf::WaveFront, path::AbstractString)
    loaded = open(path, "r") do io
        deserialize(io)
    end

    src = if loaded isa WaveFrontState
        loaded
    elseif loaded isa WaveFront
        # Backward compatibility with older state files.
        WaveFrontState(loaded)
    else
        throw(ArgumentError("State file does not contain WaveFront state"))
    end

    size(wf.field) == size(src.field) || throw(ArgumentError("State wavefront size mismatch"))

    copyto!(wf.field, src.field)
    wf.wavelength_m = src.wavelength_m
    wf.sampling_m = src.sampling_m
    wf.z_m = src.z_m
    wf.beam_diameter_m = src.beam_diameter_m
    wf.z_w0_m = src.z_w0_m
    wf.w0_m = src.w0_m
    wf.z_rayleigh_m = src.z_rayleigh_m
    wf.current_fratio = src.current_fratio
    wf.reference_surface = src.reference_surface
    wf.beam_type_old = src.beam_type_old
    wf.propagator_type = src.propagator_type
    wf.rayleigh_factor = src.rayleigh_factor
    reset_workspace!(wf.workspace)
    return wf
end
