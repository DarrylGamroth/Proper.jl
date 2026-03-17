"""Return the wavefront amplitude, `abs.(wf.field)`."""
prop_get_amplitude(wf::WaveFront) = abs.(wf.field)
