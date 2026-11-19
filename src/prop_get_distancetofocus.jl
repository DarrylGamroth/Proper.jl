"""Return the signed distance from the current plane to focus in meters."""
prop_get_distancetofocus(wf::OpticalSurfaceState) = wf.z_w0_m - wf.z_m
