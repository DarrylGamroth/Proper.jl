struct SelectPropagatorOptions{T<:AbstractFloat}
    dz::T
end

@inline SelectPropagatorOptions(wf::WaveFront{T}, dz::Real) where {T<:AbstractFloat} =
    SelectPropagatorOptions{T}(T(dz))

"""Select propagation regime based on Rayleigh-distance boundary."""
function prop_select_propagator(wf::WaveFront, opts::SelectPropagatorOptions)
    dzw = wf.z_w0_m - wf.z_m
    newz = wf.z_m + opts.dz

    beam_type_new = abs(wf.z_w0_m - newz) < wf.rayleigh_factor * wf.z_rayleigh_m ? INSIDE : OUTSIDE
    wf.propagator_type = propagator_transition(wf.beam_type_old, beam_type_new)
    wf.beam_type_old = beam_type_new
    return dzw
end

@inline function prop_select_propagator(wf::WaveFront, dz::Real)
    return prop_select_propagator(wf, SelectPropagatorOptions(wf, dz))
end
