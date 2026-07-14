using AMDGPU

AMDGPU.functional() || error("AMDGPU.functional() returned false")
AMDGPU.functional(:rocfft) || error("AMDGPU.functional(:rocfft) returned false")
AMDGPU.allowscalar(false)

include(joinpath(@__DIR__, "..", "..", "common", "prepared_latency_report.jl"))

function amdgpu_latency_field(::Type{T}, dimensions...) where {T}
    return AMDGPU.zeros(T, dimensions...)
end

function amdgpu_latency_device()
    try
        return string(AMDGPU.device())
    catch
        return "AMDGPU device"
    end
end

run_prepared_latency_report(
    :amdgpu,
    amdgpu_latency_field,
    AMDGPU.synchronize;
    device=amdgpu_latency_device(),
    backend_package_version=string(pkgversion(AMDGPU)),
)
