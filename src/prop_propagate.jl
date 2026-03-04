"""Dispatch to propagation kernel (phase-2 placeholder implementation)."""
function prop_propagate(wf::WaveFront, dz::Real, surface_name::AbstractString=""; kwargs...)
    if switch_set(:TO_PLANE; kwargs...)
        return prop_ptp(wf, dz)
    end
    return prop_ptp(wf, dz)
end
