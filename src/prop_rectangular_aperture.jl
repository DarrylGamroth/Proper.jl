function prop_rectangular_aperture(wf::WaveFront, width::Real, height::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    m = prop_rectangle(wf, width, height, xc, yc; kwargs...)
    wf.field .*= m
    return wf
end
