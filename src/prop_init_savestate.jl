"""Initialize save-state directory (creates it if needed)."""
function prop_init_savestate(path::AbstractString)
    mkpath(path)
    return path
end
