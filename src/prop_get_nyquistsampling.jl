function prop_get_nyquistsampling(wf::WaveFront, lamx::Real=0.0)
    λ = lamx == 0 ? wf.wavelength_m : float(lamx)
    return wf.current_fratio * λ / 2
end
