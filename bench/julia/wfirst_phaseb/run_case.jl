using BenchmarkTools
using JSON3
using Proper
using Statistics

@isdefined(WFIRSTPhaseBProper) || include(joinpath(@__DIR__, "..", "..", "..", "reference_models", "wfirst_phaseb_proper", "__init__.jl"))
using .WFIRSTPhaseBProper

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

has_flag(flag::String) = any(==(flag), ARGS)

function measure_samples(workload, samples::Integer)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    GC.gc()
    return run(@benchmarkable $workload() evals=1 samples=samples)
end

function trial_stats(trial::BenchmarkTools.Trial)
    est_med = median(trial)
    est_mean = mean(trial)
    est_min = minimum(trial)
    est_max = maximum(trial)
    return Dict(
        "timed" => true,
        "median_ns" => est_med.time,
        "mean_ns" => est_mean.time,
        "min_ns" => est_min.time,
        "max_ns" => est_max.time,
        "median_alloc_bytes" => est_med.memory,
        "mean_alloc_bytes" => est_mean.memory,
        "min_alloc_bytes" => est_min.memory,
        "max_alloc_bytes" => est_max.memory,
        "median_allocs" => est_med.allocs,
        "mean_allocs" => est_mean.allocs,
        "min_allocs" => est_min.allocs,
        "max_allocs" => est_max.allocs,
        "samples" => length(trial.times),
    )
end

function write_case_outputs(prefix::AbstractString, stack)
    for i in axes(stack, 1)
        prop_fits_write(prefix * "_$(i)_real.fits", Float64.(real(@view stack[i, :, :])))
        prop_fits_write(prefix * "_$(i)_imag.fits", Float64.(imag(@view stack[i, :, :])))
    end
end

function summarize_output(stack, samplings)
    mags = abs.(stack)
    return Dict(
        "dtype" => string(eltype(stack)),
        "shape" => [size(stack, 1), size(stack, 2), size(stack, 3)],
        "sampling" => Float64.(samplings),
        "peak_abs" => maximum(mags),
        "mean_abs" => sum(mags) / length(mags),
        "sum_abs2" => sum(abs2, stack),
    )
end

function main()
    case_name = String(arg_value("--case", "compact_hlc"))
    parity_only = has_flag("--parity-only")
    samples = parse(Int, String(arg_value("--samples", parity_only ? "0" : "3")))
    data_root = String(arg_value("--data-root", phaseb_default_data_root()))
    threaded = String(arg_value("--threaded", "true")) != "false"

    cases = phaseb_case_definitions()
    haskey(cases, case_name) || error("unsupported case $(case_name)")
    case = cases[case_name]
    models, _assets = prepare_phaseb_models(case; data_root=data_root)

    workload() = run_phaseb_case(models; threaded=threaded)
    timed = !parity_only && samples > 0
    workload() # warmup to exclude TTFx from timing
    trial = timed ? measure_samples(workload, samples) : nothing
    stats = timed ? trial_stats(trial) : Dict(
        "timed" => false,
        "median_ns" => nothing,
        "mean_ns" => nothing,
        "min_ns" => nothing,
        "max_ns" => nothing,
        "median_alloc_bytes" => nothing,
        "mean_alloc_bytes" => nothing,
        "min_alloc_bytes" => nothing,
        "max_alloc_bytes" => nothing,
        "median_allocs" => nothing,
        "mean_allocs" => nothing,
        "min_allocs" => nothing,
        "max_allocs" => nothing,
        "samples" => 0,
    )
    stack, samplings = workload()

    outdir = joinpath(@__DIR__, "..", "..", "reports")
    mkpath(outdir)
    prefix = joinpath(outdir, "julia_wfirst_phaseb_$(case_name)")
    write_case_outputs(prefix, stack)

    report = Dict(
        "meta" => Dict(
            "run_tag" => "wfirst_phaseb_julia_$(case_name)",
            "backend" => "cpu",
            "julia_version" => string(VERSION),
            "baseline" => "wfirst_phaseb_julia",
            "case" => case_name,
            "description" => case.description,
            "output_dim" => case.output_dim,
            "wavelengths_um" => Float64.(case.wavelengths_um),
            "available" => true,
            "data_root" => abspath(data_root),
            "threads" => Threads.nthreads(),
        ),
        "stats" => stats,
        "output" => summarize_output(stack, samplings),
        "policy" => timed ?
            "Julia CPU WFIRST Phase B reference model timing; uses BenchmarkTools with evals=1 after one untimed warmup; TTFx excluded; uses PreparedModel with shared cached assets where available" :
            "Julia CPU WFIRST Phase B reference model parity-only run; no timing samples collected",
    )

    open(prefix * ".json", "w") do io
        JSON3.write(io, report)
    end
    println(JSON3.pretty(JSON3.read(JSON3.write(report))))
end

main()
