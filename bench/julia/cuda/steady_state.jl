using BenchmarkTools
using JSON3
using Proper
include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "cuda_support.jl"))
using .BenchMetadata

const RUN_TAG = "steady_state_cuda"
const REPORT_PATH = joinpath(@__DIR__, "..", "..", "reports", "julia_cuda_steady_state.json")
const GRID_N = 512
const N_SAMPLES = 20

function workload(cuda_mod)
    wf = cuda_wavefront_begin(cuda_mod, 2.4, 0.55e-6, GRID_N; beam_diam_fraction=0.5)
    ctx = RunContext(typeof(wf.field))
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0, ctx)
    prop_end(wf)
    cuda_mod.synchronize()
    return nothing
end

function main()
    cuda_mod, reason = load_cuda_backend()
    if cuda_mod === nothing
        return write_benchmark_report(REPORT_PATH, skipped_cuda_report(RUN_TAG, reason))
    end

    # Warmup: exclude compilation and initial CUDA context setup.
    workload(cuda_mod)

    trial = run(@benchmarkable begin
        workload($cuda_mod)
        $cuda_mod.synchronize()
    end evals=1 samples=N_SAMPLES)

    report = Dict(
        "meta" => merge(cuda_report_meta(RUN_TAG, cuda_mod), Dict("grid_n" => GRID_N)),
        "policy" => "steady-state supported CUDA workload timing only; TTFx and initial CUDA context setup excluded; per-sample synchronization included",
        "stats" => trial_stats(trial),
    )

    return write_benchmark_report(REPORT_PATH, report)
end

main()
