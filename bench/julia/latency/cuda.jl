using CUDA

CUDA.functional() || error("CUDA.functional() returned false")
CUDA.allowscalar(false)

include(joinpath(@__DIR__, "..", "..", "common", "prepared_latency_report.jl"))

function cuda_latency_field(::Type{T}, dimensions...) where {T}
    return CUDA.zeros(T, dimensions...)
end

function cuda_latency_device()
    try
        return CUDA.name(CUDA.device())
    catch
        return string(CUDA.device())
    end
end

run_prepared_latency_report(
    :cuda,
    cuda_latency_field,
    CUDA.synchronize;
    device=cuda_latency_device(),
    backend_package_version=string(pkgversion(CUDA)),
)
