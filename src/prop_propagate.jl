struct PropagateOptions
    to_plane::Bool
end

@inline PropagateOptions(kwargs::Base.Iterators.Pairs) = PropagateOptions(kw_lookup_bool(kwargs, :TO_PLANE, false))

@inline function _prop_propagate!(wf::WaveFront, dz::Real, opts::PropagateOptions)
    return _prop_propagate!(wf, dz, opts, RunContext(typeof(wf.field)))
end

@inline function _prop_propagate!(wf::WaveFront, dz::Real, opts::PropagateOptions, ctx::RunContext)
    prop_select_propagator(wf, dz)

    dzv = float(dz)
    z1 = wf.z_m
    z2 = z1 + dzv

    if opts.to_plane
        if wf.propagator_type == :INSIDE__to_OUTSIDE
            wf.propagator_type = :INSIDE__to_INSIDE_
        elseif wf.propagator_type == :OUTSIDE_to_OUTSIDE
            wf.propagator_type = :OUTSIDE_to_INSIDE_
        end
    end

    if wf.propagator_type == :INSIDE__to_INSIDE_
        prop_ptp(wf, dzv, ctx)
    elseif wf.propagator_type == :INSIDE__to_OUTSIDE
        prop_ptp(wf, wf.z_w0_m - z1, ctx)
        prop_wts(wf, z2 - wf.z_w0_m, ctx)
    elseif wf.propagator_type == :OUTSIDE_to_INSIDE_
        prop_stw(wf, wf.z_w0_m - z1, ctx)
        prop_ptp(wf, z2 - wf.z_w0_m, ctx)
    elseif wf.propagator_type == :OUTSIDE_to_OUTSIDE
        prop_stw(wf, wf.z_w0_m - z1, ctx)
        prop_wts(wf, z2 - wf.z_w0_m, ctx)
    else
        throw(ArgumentError("Unknown propagator_type $(wf.propagator_type)"))
    end

    return wf
end

"""Select and execute appropriate propagation regime for distance `dz`."""
function prop_propagate(
    wf::WaveFront,
    dz::Real,
    ctx::RunContext,
    surface_name::AbstractString="";
    kwargs...,
)
    _ = surface_name
    opts = PropagateOptions(kwargs)
    return _prop_propagate!(wf, dz, opts, ctx)
end

"""Select and execute appropriate propagation regime for distance `dz`."""
function prop_propagate(
    wf::WaveFront,
    dz::Real,
    surface_name::AbstractString="";
    kwargs...,
)
    _ = surface_name
    opts = PropagateOptions(kwargs)
    return _prop_propagate!(wf, dz, opts)
end
