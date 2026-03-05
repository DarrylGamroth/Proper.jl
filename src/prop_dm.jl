"""Apply deformable mirror phase map in meters (wavefront mode by default)."""
function prop_dm(wf::WaveFront, dm_map::AbstractMatrix; mirror::Union{Nothing,Bool}=nothing, kwargs...)
    size(dm_map) == size(wf.field) || throw(ArgumentError("dm_map size must match wavefront"))
    mirror_val = mirror === nothing ? false : mirror
    # Compatibility alias: accept legacy uppercase keyword at API boundary.
    if haskey(kwargs, :MIRROR) || haskey(kwargs, :mirror)
        mirror_compat = haskey(kwargs, :MIRROR) ? compat_bool(kwargs[:MIRROR]) : compat_bool(kwargs[:mirror])
        if mirror !== nothing && mirror_compat != mirror
            throw(ArgumentError("Conflicting values for mirror and MIRROR"))
        end
        mirror_val = mirror_compat
    end
    scale = mirror_val ? 4pi / wf.wavelength_m : 2pi / wf.wavelength_m
    wf.field .*= cis.(scale .* dm_map)
    return wf
end
