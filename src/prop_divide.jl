"""Divide the current wavefront amplitude by a scalar or centered map."""
function prop_divide(wf::WaveFront, d::Number)
    wf.field ./= d
    return wf
end

"""Divide the current wavefront amplitude by a 2-D map centered at `(n/2, n/2)`."""
function prop_divide(wf::WaveFront, d::AbstractMatrix)
    size(d) == size(wf.field) || throw(ArgumentError("divisor size must match wavefront"))
    if same_backend_style(typeof(wf.field), typeof(d)) && wf.field isa StridedMatrix
        scratch = if eltype(d) <: Real
            ensure_fft_real_scratch!(wf.workspace.fft, size(d, 1), size(d, 2))
        else
            ensure_fft_scratch!(wf.workspace.fft, size(d, 1), size(d, 2))
        end
        prop_shift_center!(scratch, d; inverse=true)
        wf.field ./= scratch
    else
        wf.field ./= backend_adapt(wf.field, prop_shift_center(d; inverse=true))
    end
    return wf
end
