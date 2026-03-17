"""Return the wavefront phase, `angle.(wf.field)`."""
prop_get_phase(wf::WaveFront) = angle.(wf.field)
