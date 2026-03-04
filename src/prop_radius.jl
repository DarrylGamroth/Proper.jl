"""Return radial coordinate map (meters) at current sampling."""
function prop_radius(wf::WaveFront; kwargs...)
    ny, nx = size(wf.field)
    return radius_map(ny, nx, wf.sampling_m)
end
