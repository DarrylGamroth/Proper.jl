"""Select and execute appropriate propagation regime for distance `dz`."""
function prop_propagate(wf::WaveFront, dz::Real, surface_name::AbstractString=""; kwargs...)
    prop_select_propagator(wf, dz)

    z1 = wf.z_m
    z2 = z1 + float(dz)

    if switch_set(:TO_PLANE; kwargs...)
        if wf.propagator_type == :INSIDE__to_OUTSIDE
            wf.propagator_type = :INSIDE__to_INSIDE_
        elseif wf.propagator_type == :OUTSIDE_to_OUTSIDE
            wf.propagator_type = :OUTSIDE_to_INSIDE_
        end
    end

    if wf.propagator_type == :INSIDE__to_INSIDE_
        prop_ptp(wf, dz)
    elseif wf.propagator_type == :INSIDE__to_OUTSIDE
        prop_ptp(wf, wf.z_w0_m - z1)
        prop_wts(wf, z2 - wf.z_w0_m)
    elseif wf.propagator_type == :OUTSIDE_to_INSIDE_
        prop_stw(wf, wf.z_w0_m - z1)
        prop_ptp(wf, z2 - wf.z_w0_m)
    elseif wf.propagator_type == :OUTSIDE_to_OUTSIDE
        prop_stw(wf, wf.z_w0_m - z1)
        prop_wts(wf, z2 - wf.z_w0_m)
    else
        throw(ArgumentError("Unknown propagator_type $(wf.propagator_type)"))
    end

    return wf
end
