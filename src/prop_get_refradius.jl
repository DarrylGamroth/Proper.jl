"""Return the signed reference-sphere radius in meters."""
prop_get_refradius(wf::WaveFront) = wf.z_m - wf.z_w0_m
