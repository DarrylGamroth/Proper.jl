function prop_circular_aperture(wf::WaveFront, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    m = prop_ellipse(wf, radius, radius, xc, yc; kwargs...)
    wf.field .*= m
    return wf
end
