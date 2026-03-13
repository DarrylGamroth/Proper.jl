using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_steady_state_workload.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda_fp32"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state_fp32.json")

CUDA.allowscalar(false)
run_cuda_steady_state_report(Float32, RUN_TAG, REPORT_PATH)
