"""Spherical-to-waist propagation placeholder."""
function prop_stw(wf::WaveFront, dz::Real)
    prop_qphase(wf, -float(dz))
    wf.z_m += float(dz)
    return wf
end
