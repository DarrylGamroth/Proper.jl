using BenchmarkTools
using Proper

const AMDGPU_STEADY_GRID_N = 512
const AMDGPU_STEADY_SAMPLES = 20

function amdgpu_steady_state_prescription(λm, n; kwargs...)
    wf = prop_begin(2.4, λm, n; beam_diam_fraction=0.5)
    prop_circular_aperture(wf, 0.6)
    prop_lens(wf, 20.0)
    prop_propagate(wf, 20.0)
    return wf
end

function amdgpu_steady_state_prepared(::Type{T}, grid_n::Integer=AMDGPU_STEADY_GRID_N) where {T<:AbstractFloat}
    wf = amdgpu_wavefront_begin(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=T(0.5))
    ctx = RunContext(wf)
    return prepare_model(:steady_state_amdgpu, amdgpu_steady_state_prescription, T(0.55), grid_n; context=ctx, pool_size=1)
end

function amdgpu_steady_state_workload(prepared::Union{PreparedPrescription,PreparedModel})
    prop_run(prepared)
    amdgpu_sync()
    return nothing
end

function run_amdgpu_steady_state_report(
    ::Type{T},
    run_tag::String,
    report_path::AbstractString;
    grid_n::Integer=AMDGPU_STEADY_GRID_N,
    samples::Integer=AMDGPU_STEADY_SAMPLES,
) where {T<:AbstractFloat}
    prepared = amdgpu_steady_state_prepared(T, grid_n)
    amdgpu_steady_state_workload(prepared)

    trial = run(@benchmarkable begin
        amdgpu_steady_state_workload($prepared)
        amdgpu_sync()
    end evals=1 samples=samples)

    report = Dict(
        "meta" => merge(amdgpu_report_meta(run_tag; device=amdgpu_device_label()), Dict("grid_n" => grid_n, "precision" => string(T))),
        "policy" => "steady-state supported AMDGPU workload timing only; TTFx and initial GPU context setup excluded; per-sample synchronization included",
        "stats" => trial_stats(trial),
    )

    return write_benchmark_report(report_path, report)
end
