function prop_multiply(wf::WaveFront, m::Number)
    wf.field .*= m
    return wf
end

function prop_multiply(wf::WaveFront, m::AbstractMatrix)
    size(m) == size(wf.field) || throw(ArgumentError("multiplier size must match wavefront"))
    wf.field .*= m
    return wf
end
