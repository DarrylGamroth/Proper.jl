"""Apply deformable mirror phase map in meters (wavefront mode by default)."""
function prop_dm(wf::WaveFront, dm_map::AbstractMatrix; MIRROR::Union{Bool,Integer}=false)
    size(dm_map) == size(wf.field) || throw(ArgumentError("dm_map size must match wavefront"))
    scale = compat_bool(MIRROR) ? 4pi / wf.wavelength_m : 2pi / wf.wavelength_m
    wf.field .*= cis.(scale .* dm_map)
    return wf
end
