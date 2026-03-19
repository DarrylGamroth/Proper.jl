using BenchmarkTools
using Proper
include(joinpath(@__DIR__, "prepared_execution_workloads.jl"))

const CUDA_STEADY_GRID_N = 512
const CUDA_STEADY_SAMPLES = 20

function cuda_steady_state_prepared(::Type{T}, grid_n::Integer=CUDA_STEADY_GRID_N) where {T<:AbstractFloat}
    ctx = backend_prepared_context(cuda_wavefront_begin, T, grid_n)
    return prepare_steady_state_model(T, grid_n, ctx; name=:steady_state_cuda)
end

function cuda_steady_state_workload(prepared::Union{PreparedPrescription,PreparedModel})
    return run_prepared_workload(prepared, cuda_sync)
end

function _run_cuda_steady_state_report(::Type{T}, run_tag::String, report_path::AbstractString; grid_n::Integer=CUDA_STEADY_GRID_N, samples::Integer=CUDA_STEADY_SAMPLES) where {T<:AbstractFloat}
    prepared = cuda_steady_state_prepared(T, grid_n)
    trial = benchmark_prepared_trial(() -> cuda_steady_state_workload(prepared); samples=samples)

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
