"""Divide the current wavefront amplitude by a scalar or centered map."""
function prop_divide(wf::WaveFront, d::Number)
    wf.field ./= d
    return wf
end

"""Divide the current wavefront amplitude by a 2-D map centered at `(n/2, n/2)`."""
function prop_divide(wf::WaveFront, d::AbstractMatrix)
    size(d) == size(wf.field) || throw(ArgumentError("divisor size must match wavefront"))
    wf.field ./= shift_center_for_wavefront!(wf, d; inverse=true)
    return wf
end
