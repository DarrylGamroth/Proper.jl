using Proper

@isdefined(WFIRSTPhaseBProper) || include(joinpath(@__DIR__, "..", "..", "..", "reference_models", "wfirst_phaseb_proper", "__init__.jl"))
using .WFIRSTPhaseBProper

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "benchmark_helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "suites", "performance", "reporting.jl"))

using .BenchMetadata
using .PerformanceBenchmarkHelpers
using .BenchmarkReporting

function bench_wfirst_case(case_name::AbstractString, data_root::AbstractString, samples::Integer, threaded::Bool, skip_missing::Bool)
    cases = phaseb_case_definitions()
    haskey(cases, case_name) || throw(ArgumentError("unsupported WFIRST Phase B case: $(case_name)"))
    case = cases[case_name]

    try
        models, _assets = prepare_phaseb_models(case; data_root=data_root)
        workload() = begin
            run_phaseb_case(models; threaded=threaded)
            nothing
        end
        stats = trial_stats(measure_samples(workload, samples))
        return Dict(
            "case" => String(case_name),
            "available" => true,
            "description" => case.description,
            "output_dim" => case.output_dim,
            "wavelengths_um" => Float64.(case.wavelengths_um),
            "threaded" => threaded,
            "stats" => stats,
        )
    catch err
        skip_missing || rethrow()
        return Dict(
            "case" => String(case_name),
            "available" => false,
            "description" => case.description,
            "output_dim" => case.output_dim,
            "wavelengths_um" => Float64.(case.wavelengths_um),
            "threaded" => threaded,
            "stats" => untimed_stats(),
            "error" => sprint(showerror, err),
        )
    end
end

function main()
    samples = parse(Int, String(arg_value("--samples", "3")))
    data_root = String(arg_value("--data-root", phaseb_default_data_root()))
    threaded = String(arg_value("--threaded", "true")) != "false"
    skip_missing = String(arg_value("--skip-missing", "true")) != "false"
    case_names = parse_string_list(arg_value("--cases", nothing), DEFAULT_WFIRST_CASES)

    cases = [bench_wfirst_case(case_name, data_root, samples, threaded, skip_missing) for case_name in case_names]
    report = Dict(
        "meta" => merge(
            benchmark_metadata(run_tag="wfirst_prepared_models", backend=:cpu, baseline="wfirst_phaseb_julia"),
            Dict("benchmark" => "wfirst_prepared_models", "data_root" => abspath(data_root)),
        ),
        "policy" => "Steady-state CPU WFIRST/Roman prepared-model benchmarks via BenchmarkTools with evals=1 after warmup; TTFx excluded. Data assets are read from the external cache path and are not bundled into the package.",
        "cases" => cases,
    )
    out = joinpath(@__DIR__, "..", "..", "reports", "wfirst_prepared_models.json")
    write_json_report(out, report)
end

main()
