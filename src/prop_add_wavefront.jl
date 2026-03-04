"""Add complex wavefront samples."""
function prop_add_wavefront(wf::WaveFront, other::AbstractMatrix)
    size(other) == size(wf.field) || throw(ArgumentError("wavefront size must match"))
    wf.field .+= other
    return wf
end
