using BenchmarkTools
using JSON3
using Proper
using CUDA

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_wavefront_kernel_cases.jl"))
using .BenchMetadata

const RUN_TAG = "cuda_isolated_wavefront_kernel"
const GRID_N = 512
const WARMUP_ITERS = 40
const N_SAMPLES = 60

CUDA.allowscalar(false)

function parse_precision(label::AbstractString)
    label == "fp64" && return Float64
    label == "fp32" && return Float32
    throw(ArgumentError("precision must be fp64 or fp32"))
end

precision_label(::Type{Float64}) = "fp64"
precision_label(::Type{Float32}) = "fp32"

function report_path(kernel::AbstractString, precision::AbstractString)
    return joinpath(@__DIR__, "..", "..", "reports", "cuda_isolated_$(kernel)_$(precision).json")
end

function benchmark_isolated_wavefront_kernel(kernel::AbstractString, ::Type{T}) where {T<:AbstractFloat}
    Symbol(kernel) in Symbol.(CUDA_WAVEFRONT_KERNEL_ORDER) ||
        throw(ArgumentError("unsupported kernel '$kernel'"))

    cases = build_cuda_wavefront_kernel_cases(T, GRID_N)
    case = getproperty(cases, Symbol(kernel))

    warmup_cuda_benchmark_case(case, WARMUP_ITERS)
    host_trial = run_cuda_benchmark_case(case, N_SAMPLES)
    timing_stats = collect_cuda_host_device_samples(case, N_SAMPLES)

    return Dict(
        "meta" => merge(
            cuda_report_meta(RUN_TAG; device=cuda_device_label()),
            Dict(
                "kernel" => kernel,
                "precision" => precision_label(T),
                "grid_n" => GRID_N,
                "warmup_iters" => WARMUP_ITERS,
            ),
        ),
        "policy" => "isolated one-kernel-per-process CUDA microbenchmark; host wall time uses BenchmarkTools with per-sample state restore, device time uses CUDA.@elapsed with explicit pre-sample synchronization; TTFx excluded",
        "host" => trial_stats(host_trial),
        "timing" => timing_stats,
    )
end

function main(args)
    length(args) == 2 || throw(ArgumentError("usage: isolated_wavefront_kernel.jl <kernel> <fp64|fp32>"))
    kernel = args[1]
    precision = args[2]
    T = parse_precision(precision)
    return write_benchmark_report(report_path(kernel, precision), benchmark_isolated_wavefront_kernel(kernel, T))
end

main(ARGS)
