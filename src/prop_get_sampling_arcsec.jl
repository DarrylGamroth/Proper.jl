const _RAD_TO_ARCSEC = 206264.80624709636
prop_get_sampling_arcsec(wf::WaveFront) = prop_get_sampling_radians(wf) * _RAD_TO_ARCSEC
