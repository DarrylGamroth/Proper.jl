"""Placeholder DM fit returning input map and zero residual."""
function prop_fit_dm(dm_map::AbstractMatrix; kwargs...)
    return copy(dm_map), zero(eltype(dm_map))
end
