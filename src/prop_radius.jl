"""Return radial coordinate map at current sampling; `NORM` normalizes by beam radius."""
function prop_radius(wf::WaveFront; kwargs...)
    ny, nx = size(wf.field)
    r = radius_map(ny, nx, wf.sampling_m)
    if switch_set(:NORM; kwargs...)
        return r ./ prop_get_beamradius(wf)
    end
    return r
end
