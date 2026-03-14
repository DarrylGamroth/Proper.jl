using BenchmarkTools
using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_wavefront_kernel_cases.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda_precision_split"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "cuda_precision_split.json")
const GRID_N = 512
const N_SAMPLES = 20

CUDA.allowscalar(false)

function benchmark_precision_split()
    stats64 = benchmark_cuda_wavefront_kernel_stats(Float64; grid_n=GRID_N, samples=N_SAMPLES)
    stats32 = benchmark_cuda_wavefront_kernel_stats(Float32; grid_n=GRID_N, samples=N_SAMPLES)

    report = Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG; device=cuda_device_label()), Dict("grid_n" => GRID_N)),
        "policy" => "CUDA precision-split timing for steady-state propagation-heavy workload and selected kernels with per-sample wavefront state restore; TTFx excluded; per-sample synchronization included",
        "kernels" => Dict(
            "prop_qphase" => Dict("fp64" => stats64["prop_qphase"], "fp32" => stats32["prop_qphase"]),
            "prop_ptp" => Dict("fp64" => stats64["prop_ptp"], "fp32" => stats32["prop_ptp"]),
            "prop_wts" => Dict("fp64" => stats64["prop_wts"], "fp32" => stats32["prop_wts"]),
            "prop_stw" => Dict("fp64" => stats64["prop_stw"], "fp32" => stats32["prop_stw"]),
            "prop_circular_aperture" => Dict("fp64" => stats64["prop_circular_aperture"], "fp32" => stats32["prop_circular_aperture"]),
            "prop_end_mutating" => Dict("fp64" => stats64["prop_end_mutating"], "fp32" => stats32["prop_end_mutating"]),
        ),
    )

    return write_benchmark_report(REPORT_PATH, report)
end

benchmark_precision_split()
