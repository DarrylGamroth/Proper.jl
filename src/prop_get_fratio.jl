prop_get_fratio(wf::WaveFront) = wf.beam_diameter_m == 0 ? Inf : abs(wf.z_m) / wf.beam_diameter_m
