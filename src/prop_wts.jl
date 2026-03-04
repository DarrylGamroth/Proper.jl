"""Waist-to-spherical propagation placeholder."""
function prop_wts(wf::WaveFront, dz::Real)
    prop_qphase(wf, float(dz))
    wf.z_m += float(dz)
    return wf
end
