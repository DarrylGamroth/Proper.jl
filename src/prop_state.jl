using Serialization

"""Load serialized state and copy fields into `wf` (Python bugfix behavior)."""
function prop_state(wf::WaveFront, path::AbstractString)
    loaded = open(path, "r") do io
        deserialize(io)
    end
    loaded isa WaveFront || throw(ArgumentError("State file does not contain WaveFront"))

    wf.field .= loaded.field
    wf.wavelength_m = loaded.wavelength_m
    wf.sampling_m = loaded.sampling_m
    wf.z_m = loaded.z_m
    wf.beam_diameter_m = loaded.beam_diameter_m
    return wf
end
