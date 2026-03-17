"""Return the current focal-plane sampling in radians per pixel."""
function prop_get_sampling_radians(wf::WaveFront)
    fl = prop_get_fratio(wf) * wf.beam_diameter_m
    return prop_get_sampling(wf) / fl
end
