"""Multiply the current wavefront amplitude by a scalar or centered map."""
function prop_multiply(wf::WaveFront, m::Number)
    wf.field .*= m
    return wf
end

"""Multiply the current wavefront amplitude by a 2-D map centered at `(n/2, n/2)`."""
function prop_multiply(wf::WaveFront, m::AbstractMatrix)
    size(m) == size(wf.field) || throw(ArgumentError("multiplier size must match wavefront"))
    if same_backend_style(typeof(wf.field), typeof(m)) && wf.field isa StridedMatrix
        scratch = if eltype(m) <: Real
            ensure_fft_real_scratch!(wf.workspace.fft, size(m, 1), size(m, 2))
        else
            ensure_fft_scratch!(wf.workspace.fft, size(m, 1), size(m, 2))
        end
        prop_shift_center!(scratch, m; inverse=true)
        wf.field .*= scratch
    else
        wf.field .*= backend_adapt(wf.field, prop_shift_center(m; inverse=true))
    end
    return wf
end
