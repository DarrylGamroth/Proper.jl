const _RAD_TO_ARCSEC = 206264.80624709636

"""Return the current focal-plane sampling in arcseconds per pixel."""
function prop_get_sampling_arcsec(wf::WaveFront)
    fl = prop_get_fratio(wf) * wf.beam_diameter_m
    return prop_get_sampling(wf) * _RAD_TO_ARCSEC / fl
end
