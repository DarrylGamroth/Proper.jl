"""Return the wavefront grid size in pixels."""
prop_get_gridsize(wf::WaveFront) = size(wf.field, 1)
prop_get_gridsize(wf::WaveFrontGeometry) = wf.gridsize
