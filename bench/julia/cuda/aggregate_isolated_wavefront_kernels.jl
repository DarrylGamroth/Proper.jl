using JSON3

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_wavefront_kernel_cases.jl"))
using .BenchMetadata

const RUN_TAG = "cuda_isolated_wavefront_kernels"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_isolated_wavefront_kernels.json")

function kernel_report_path(kernel::AbstractString, precision::AbstractString)
    return joinpath(@__DIR__, "..", "..", "reports", "cuda_isolated_$(kernel)_$(precision).json")
end

function load_report(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing isolated CUDA kernel report: $path"))
    return JSON3.read(read(path, String))
end

function getkey(obj, key)
    return haskey(obj, key) ? obj[key] : nothing
end

function aggregate_reports()
    kernels = Dict{String,Any}()
    first_device = nothing

    for kernel in CUDA_WAVEFRONT_KERNEL_ORDER
        payload = Dict{String,Any}()
        for precision in ("fp64", "fp32")
            report = load_report(kernel_report_path(kernel, precision))
            first_device === nothing && (first_device = getkey(getkey(report, "meta"), "device"))
            payload[precision] = Dict(
                "host" => getkey(report, "host"),
                "timing" => getkey(report, "timing"),
            )
        end
        kernels[kernel] = payload
    end

    meta = first_device === nothing ? cuda_report_meta(RUN_TAG) : cuda_report_meta(RUN_TAG; device=String(first_device))
    return Dict(
        "meta" => meta,
        "policy" => "isolated one-kernel-per-process CUDA wavefront microbenchmarks with longer warmup and both host wall and device timing",
        "kernels" => kernels,
    )
end

write_benchmark_report(REPORT_PATH, aggregate_reports())
