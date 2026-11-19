"""Return the current plane sampling in meters per pixel."""
prop_get_sampling(wf::OpticalSurfaceState) = wf.sampling_m
