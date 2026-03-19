using BenchmarkTools
using Proper
include(joinpath(@__DIR__, "prepared_execution_workloads.jl"))

const AMDGPU_STEADY_GRID_N = 512
const AMDGPU_STEADY_SAMPLES = 20

function amdgpu_steady_state_prepared(::Type{T}, grid_n::Integer=AMDGPU_STEADY_GRID_N) where {T<:AbstractFloat}
    ctx = backend_prepared_context(amdgpu_wavefront_begin, T, grid_n)
    return prepare_steady_state_model(T, grid_n, ctx; name=:steady_state_amdgpu)
end

function amdgpu_steady_state_workload(prepared::Union{PreparedPrescription,PreparedModel})
    return run_prepared_workload(prepared, amdgpu_sync)
end

function run_amdgpu_steady_state_report(
    ::Type{T},
    run_tag::String,
    report_path::AbstractString;
    grid_n::Integer=AMDGPU_STEADY_GRID_N,
    samples::Integer=AMDGPU_STEADY_SAMPLES,
    alias_specs=(),
) where {T<:AbstractFloat}
    prepared = amdgpu_steady_state_prepared(T, grid_n)
    trial = benchmark_prepared_trial(() -> amdgpu_steady_state_workload(prepared); samples=samples)

    report = Dict(
        "meta" => merge(amdgpu_report_meta(run_tag; device=amdgpu_device_label()), Dict("grid_n" => grid_n, "precision" => string(T))),
        "policy" => "steady-state supported AMDGPU workload timing only; TTFx and initial GPU context setup excluded; per-sample synchronization included",
        "stats" => trial_stats(trial),
    )

    write_benchmark_report(report_path, report)
    for alias in alias_specs
        path, alias_run_tag = alias
        aliased = copy(report)
        meta = copy(report["meta"])
        meta["run_tag"] = alias_run_tag
        aliased["meta"] = meta
        write_benchmark_report(path, aliased)
    end
    return report
end
