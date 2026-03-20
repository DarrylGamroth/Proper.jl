"""
    prop_end_savestate(wf, path; kwargs...)

Save state then finalize wavefront output.

# Notes
- State serialization is currently host-only. `wf.field` must use a CPU-backed
  array backend.
"""
function prop_end_savestate(wf::WaveFront, path::AbstractString; kwargs...)
    prop_savestate(wf, path)
    return prop_end(wf; kwargs...)
end
