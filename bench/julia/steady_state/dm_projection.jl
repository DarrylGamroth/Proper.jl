using Proper

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "benchmark_helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "reporting.jl"))

using .BenchMetadata
using .PerformanceBenchmarkHelpers
using .BenchmarkReporting

function bench_dm_case(active_count::Integer, grid_n::Integer, samples::Integer)
    cmd = active_dm_command(active_count)
    count_active_actuators(cmd.values) == cmd.active_count || error("active actuator count mismatch")

    wf_map = prop_begin(1.0, 0.55e-6, grid_n)
    map_workload() = begin
        prop_dm(wf_map, cmd.values, 0.0, 0.0, 0.0; N_ACT_ACROSS_PUPIL=cmd.grid_side, NO_APPLY=true)
        nothing
    end

    wf_apply = prop_begin(1.0, 0.55e-6, grid_n)
    unit_field = one(eltype(wf_apply.field))
    apply_workload() = begin
        fill!(wf_apply.field, unit_field)
        prop_dm(wf_apply, cmd.values, 0.0, 0.0, 0.0; N_ACT_ACROSS_PUPIL=cmd.grid_side)
        nothing
    end

    return Dict(
        "label" => cmd.label,
        "active_count" => cmd.active_count,
        "storage_grid_side" => cmd.grid_side,
        "storage_grid" => "$(cmd.grid_side)x$(cmd.grid_side)",
        "wavefront_grid" => grid_n,
        "no_apply_map" => trial_stats(measure_samples(map_workload, samples)),
        "applied_wavefront" => trial_stats(measure_samples(apply_workload, samples)),
    )
end

function main()
    samples = parse(Int, String(arg_value("--samples", "5")))
    grid_n = parse(Int, String(arg_value("--grid", "256")))
    counts = parse_int_list(arg_value("--active-counts", nothing), DEFAULT_DM_ACTIVE_COUNTS)

    cases = [bench_dm_case(active_count, grid_n, samples) for active_count in counts]
    report = Dict(
        "meta" => merge(
            benchmark_metadata(run_tag="dm_projection", backend=:cpu),
            Dict("benchmark" => "dm_projection"),
        ),
        "policy" => "Steady-state CPU DM projection benchmarks via BenchmarkTools with evals=1 after warmup; TTFx excluded. Active actuator counts are embedded in the smallest square storage grid with centered deterministic active masks.",
        "cases" => cases,
    )
    out = joinpath(@__DIR__, "..", "..", "reports", "dm_projection.json")
    write_json_report(out, report)
end

main()
