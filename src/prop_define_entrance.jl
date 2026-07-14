@inline _entrance_power(::CPUBackend, field::AbstractArray) = sum(abs2, field)
@inline _entrance_power(::UnknownBackend, field::AbstractArray) = sum(abs2, field)

@inline function _entrance_power(::Union{CUDABackend,AMDGPUBackend}, field::AbstractArray)
    T = real(eltype(field))
    neutral = zero(T)
    return AK.mapreduce(abs2, +, field; init=neutral, neutral=neutral)
end

"""
    prop_define_entrance(wf)

Normalize a finite, nonzero entrance pupil to unit total intensity.

The wavefront backend is preserved. A zero or nonfinite total intensity is
rejected because no finite normalization exists for that input.
"""
function prop_define_entrance(wf::WaveFront)
    p = _entrance_power(backend_style(typeof(wf.field)), wf.field)
    isfinite(p) && p > zero(p) || throw(DomainError(
        p,
        "prop_define_entrance requires finite positive total intensity",
    ))
    wf.field .*= inv(sqrt(p))
    return wf
end
