"""Hex wavefront placeholder returns circular aperture equivalent."""
function prop_hex_wavefront(wf::WaveFront, radius::Real)
    return prop_circular_aperture(wf, radius)
end
