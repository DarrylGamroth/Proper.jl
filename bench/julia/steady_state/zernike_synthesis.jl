using Proper

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "benchmark_helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "reporting.jl"))

using .BenchMetadata
using .PerformanceBenchmarkHelpers
using .BenchmarkReporting

function bench_zernike_synthesis_case(grid_n::Integer, nterms::Integer, samples::Integer)
    wf = prop_begin(1.0, 0.55e-6, grid_n)
    nums = collect(1:Int(nterms))
    coeffs = zernike_coefficients(nterms)
    workload() = begin
        prop_zernikes(wf, nums, coeffs; no_apply=true)
        nothing
    end
    return Dict(
        "grid" => Int(grid_n),
        "nterms" => Int(nterms),
        "variant" => "weighted_sum_no_apply",
        "stats" => trial_stats(measure_samples(workload, samples)),
    )
end

function main()
    samples = parse(Int, String(arg_value("--samples", "5")))
    grids = parse_int_list(arg_value("--grids", nothing), DEFAULT_ZERNIKE_GRIDS)
    nterms_list = parse_int_list(arg_value("--nterms", nothing), DEFAULT_ZERNIKE_NTERMS)

    cases = [
        bench_zernike_synthesis_case(grid_n, nterms, samples)
        for grid_n in grids for nterms in nterms_list
    ]
    report = Dict(
        "meta" => merge(
            benchmark_metadata(run_tag="zernike_synthesis", backend=:cpu),
            Dict("benchmark" => "zernike_synthesis"),
        ),
        "policy" => "Steady-state CPU prop_zernikes(...; no_apply=true) benchmarks via BenchmarkTools with evals=1 after warmup; TTFx excluded.",
        "cases" => cases,
    )
    out = joinpath(@__DIR__, "..", "..", "reports", "zernike_synthesis.json")
    write_json_report(out, report)
end

main()
