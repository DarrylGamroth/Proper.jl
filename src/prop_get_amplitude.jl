"""Return a centered wavefront-amplitude array."""
prop_get_amplitude(wf::WaveFront) = prop_shift_center(abs.(wf.field))
