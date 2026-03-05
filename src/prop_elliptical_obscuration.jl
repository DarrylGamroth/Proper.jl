function prop_elliptical_obscuration(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    m = prop_ellipse(wf, xradius, yradius, xc, yc; kwargs...)
    wf.field .*= prop_shift_center(1 .- m)
    return wf
end
