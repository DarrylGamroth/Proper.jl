function prop_get_beamradius(wf::WaveFront)
    return wf.w0_m * sqrt(1 + (wf.wavelength_m * (wf.z_m - wf.z_w0_m) / (pi * wf.w0_m^2))^2)
end
