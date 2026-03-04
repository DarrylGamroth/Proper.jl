"""Save state then finalize wavefront output."""
function prop_end_savestate(wf::WaveFront, path::AbstractString; kwargs...)
    prop_savestate(wf, path)
    return prop_end(wf; kwargs...)
end
