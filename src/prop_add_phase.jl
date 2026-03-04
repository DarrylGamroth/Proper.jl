"""Apply a phase screen in radians."""
function prop_add_phase(wf::WaveFront, phase::AbstractMatrix)
    size(phase) == size(wf.field) || throw(ArgumentError("phase size must match wavefront"))
    @inbounds wf.field .*= cis.(phase)
    return wf
end
