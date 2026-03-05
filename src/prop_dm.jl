"""Apply deformable mirror phase map in meters (wavefront mode by default)."""
function prop_dm(wf::WaveFront, dm_map::AbstractMatrix; mirror::Bool=false)
    size(dm_map) == size(wf.field) || throw(ArgumentError("dm_map size must match wavefront"))
    scale = mirror ? 4pi / wf.wavelength_m : 2pi / wf.wavelength_m
    wf.field .*= cis.(scale .* dm_map)
    return wf
end
