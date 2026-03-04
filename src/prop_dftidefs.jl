"""Return placeholder DFTI constant dictionary for compatibility."""
function prop_dftidefs()
    return Dict(
        :DFTI_FORWARD_SCALE => 1.0,
        :DFTI_BACKWARD_SCALE => 1.0,
    )
end
