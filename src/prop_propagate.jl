struct PropagateOptions
    to_plane::Bool
end

@inline PropagateOptions(kwargs::Base.Iterators.Pairs) = PropagateOptions(kw_lookup_bool(kwargs, :TO_PLANE, false))
@inline PropagateOptions(::Base.Pairs{Symbol,Union{},Nothing,@NamedTuple{}}) = PropagateOptions(false)

@inline function _propagate_by_transition!(
    wf::WaveFront,
    z1::Real,
    z2::Real,
    dzv::Real,
    ctx::RunContext,
    pt::PropagatorType,
)
    if pt === INSIDE_TO_INSIDE
        return prop_ptp(wf, dzv, ctx)
    elseif pt === INSIDE_TO_OUTSIDE
        prop_ptp(wf, wf.z_w0_m - z1, ctx)
        return prop_wts(wf, z2 - wf.z_w0_m, ctx)
    elseif pt === OUTSIDE_TO_INSIDE
        prop_stw(wf, wf.z_w0_m - z1, ctx)
        return prop_ptp(wf, z2 - wf.z_w0_m, ctx)
    end
    prop_stw(wf, wf.z_w0_m - z1, ctx)
    return prop_wts(wf, z2 - wf.z_w0_m, ctx)
end

@inline function _prop_propagate!(wf::WaveFront, dz::Real, opts::PropagateOptions)
    return _prop_propagate!(wf, dz, opts, RunContext(wf))
end

@inline function _prop_propagate!(wf::WaveFront, dz::Real, opts::PropagateOptions, ctx::RunContext)
    prop_select_propagator(wf, dz)

    dzv = float(dz)
    z1 = wf.z_m
    z2 = z1 + dzv

    if opts.to_plane
        wf.propagator_type = to_plane_propagator(wf.propagator_type)
    end

    _propagate_by_transition!(wf, z1, z2, dzv, ctx, wf.propagator_type)

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
