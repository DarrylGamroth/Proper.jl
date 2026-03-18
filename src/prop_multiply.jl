"""Multiply the current wavefront amplitude by a scalar or centered map."""
function prop_multiply(wf::WaveFront, m::Number)
    wf.field .*= m
    return wf
end

"""Multiply the current wavefront amplitude by a 2-D map centered at `(n/2, n/2)`."""
function prop_multiply(wf::WaveFront, m::AbstractMatrix)
    size(m) == size(wf.field) || throw(ArgumentError("multiplier size must match wavefront"))
    wf.field .*= shift_center_for_wavefront!(wf, m; inverse=true)
    return wf
end
