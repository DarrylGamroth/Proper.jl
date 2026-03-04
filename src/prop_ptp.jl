"""Planar-to-planar propagation placeholder: updates z only."""
function prop_ptp(wf::WaveFront, dz::Real)
    wf.z_m += float(dz)
    return wf
end
