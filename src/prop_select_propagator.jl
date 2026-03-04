"""Return propagation mode identifier for compatibility plumbing."""
function prop_select_propagator(wf::WaveFront, dz::Real)
    return :INSIDE__to_INSIDE_
end
