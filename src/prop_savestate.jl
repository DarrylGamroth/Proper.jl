using Serialization

"""Serialize wavefront state to disk."""
function prop_savestate(wf::WaveFront, path::AbstractString)
    open(path, "w") do io
        serialize(io, wf)
    end
    return path
end
