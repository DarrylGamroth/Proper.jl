function prop_divide(wf::WaveFront, d::Number)
    wf.field ./= d
    return wf
end

function prop_divide(wf::WaveFront, d::AbstractMatrix)
    size(d) == size(wf.field) || throw(ArgumentError("divisor size must match wavefront"))
    wf.field ./= backend_adapt(wf.field, prop_shift_center(d))
    return wf
end
