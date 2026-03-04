"""Apply an 8th-order radial apodization mask placeholder."""
function prop_8th_order_mask(wf::WaveFront, radius::Real)
    r = prop_radius(wf)
    rr = r ./ float(radius)
    mask = @. exp(-(rr^8))
    wf.field .*= mask
    return wf
end
