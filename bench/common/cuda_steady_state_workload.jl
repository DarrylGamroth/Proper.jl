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

function _run_cuda_steady_state_report(::Type{T}, run_tag::String, report_path::AbstractString; grid_n::Integer=CUDA_STEADY_GRID_N, samples::Integer=CUDA_STEADY_SAMPLES) where {T<:AbstractFloat}
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

function aliased_cuda_report(report, run_tag::AbstractString)
    aliased = copy(report)
    meta = copy(report["meta"])
    meta["run_tag"] = run_tag
    aliased["meta"] = meta
    return aliased
end

function write_cuda_report_aliases(report, alias_specs::Pair{<:AbstractString,<:AbstractString}...)
    for alias in alias_specs
        path, run_tag = alias
        write_benchmark_report(path, aliased_cuda_report(report, run_tag))
    end
    return report
end

function run_cuda_steady_state_report(
    ::Type{T},
    run_tag::String,
    report_path::AbstractString;
    grid_n::Integer=CUDA_STEADY_GRID_N,
    samples::Integer=CUDA_STEADY_SAMPLES,
    alias_specs=(),
) where {T<:AbstractFloat}
    report = _run_cuda_steady_state_report(T, run_tag, report_path; grid_n=grid_n, samples=samples)
    isempty(alias_specs) || write_cuda_report_aliases(report, alias_specs...)
    return report
end
