using BenchmarkTools
using JSON3
using Proper
using CUDA
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state.json")
const GRID_N = 512
const N_SAMPLES = 20

CUDA.allowscalar(false)

function workload()
    wf = cuda_wavefront_begin(2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx = RunContext(wf)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0, ctx)
    prop_end(wf)
    cuda_sync()
    return nothing
end

function main()
    # Warmup: exclude compilation and initial CUDA context setup.
    workload()

    trial = run(@benchmarkable begin
        workload()
        cuda_sync()
    end evals=1 samples=N_SAMPLES)

    report = Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG; device=cuda_device_label()), Dict("grid_n" => GRID_N)),
        "policy" => "steady-state supported CUDA workload timing only; TTFx and initial CUDA context setup excluded; per-sample synchronization included",
        "stats" => trial_stats(trial),
    )

    return write_benchmark_report(REPORT_PATH, report)
end

main()
