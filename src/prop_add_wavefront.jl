"""Add a scalar complex amplitude to every internal wavefront sample."""
function prop_add_wavefront(wf::WaveFront, other::Number)
    wf.field .+= other
    return wf
end

"""Add a centered complex wavefront map to the internal wavefront field."""
function prop_add_wavefront(wf::WaveFront, other::AbstractMatrix)
    size(other) == size(wf.field) || throw(ArgumentError("wavefront size must match"))
    centered_to_internal = prop_shift_center(other; inverse=true)
    wf.field .+= backend_adapt(wf.field, centered_to_internal)
    return wf
end
