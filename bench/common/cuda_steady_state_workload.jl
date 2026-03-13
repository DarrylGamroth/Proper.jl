using BenchmarkTools
using Proper

const CUDA_STEADY_GRID_N = 512
const CUDA_STEADY_SAMPLES = 20

function cuda_steady_state_workload(::Type{T}, grid_n::Integer=CUDA_STEADY_GRID_N) where {T<:AbstractFloat}
    wf = cuda_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=T(0.5))
    ctx = RunContext(wf)
    prop_circular_aperture(wf, T(0.6))
    prop_lens(wf, T(20.0), ctx)
    prop_propagate(wf, T(20.0), ctx)
    prop_end(wf)
    cuda_sync()
    return nothing
end

function run_cuda_steady_state_report(::Type{T}, run_tag::String, report_path::AbstractString; grid_n::Integer=CUDA_STEADY_GRID_N, samples::Integer=CUDA_STEADY_SAMPLES) where {T<:AbstractFloat}
    cuda_steady_state_workload(T, grid_n)

    trial = run(@benchmarkable begin
        cuda_steady_state_workload($T, $grid_n)
        cuda_sync()
    end evals=1 samples=samples)

    report = Dict(
        "meta" => merge(cuda_report_meta(run_tag; device=cuda_device_label()), Dict("grid_n" => grid_n, "precision" => string(T))),
        "policy" => "steady-state supported CUDA workload timing only; TTFx and initial CUDA context setup excluded; per-sample synchronization included",
        "stats" => trial_stats(trial),
    )

    return write_benchmark_report(report_path, report)
end
