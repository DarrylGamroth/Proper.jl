"""Return a centered wavefront-phase array in radians."""
prop_get_phase(wf::WaveFront) = prop_shift_center(angle.(wf.field))
