function prop_circular_obscuration(wf::WaveFront, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    m = prop_ellipse(wf, radius, radius, xc, yc; kwargs...)
    wf.field .*= (1 .- m)
    return wf
end
