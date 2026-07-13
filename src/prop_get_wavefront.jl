"""Return a centered copy of the complex wavefront field."""
prop_get_wavefront(wf::WaveFront) = prop_shift_center(wf.field)
