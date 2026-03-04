function prop_elliptical_aperture(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    m = prop_ellipse(wf, xradius, yradius, xc, yc; kwargs...)
    wf.field .*= m
    return wf
end
