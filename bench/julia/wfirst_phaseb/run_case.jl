using FITSIO
using JSON3
using Proper
using Statistics
using Proper.WFIRSTPhaseBProper

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

function time_samples_ns(workload, samples::Integer)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    times = Vector{Int64}(undef, samples)
    GC.gc()
    for i in 1:samples
        t0 = time_ns()
        workload()
        times[i] = time_ns() - t0
    end
    return times
end

function write_case_outputs(prefix::AbstractString, stack)
    for i in axes(stack, 1)
        FITS(prefix * "_$(i)_real.fits", "w") do f
            write(f, Float64.(real(@view stack[i, :, :])))
        end
        FITS(prefix * "_$(i)_imag.fits", "w") do f
            write(f, Float64.(imag(@view stack[i, :, :])))
        end
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
    samples = parse(Int, String(arg_value("--samples", "3")))
    data_root = String(arg_value("--data-root", phaseb_default_data_root()))
    threaded = String(arg_value("--threaded", "true")) != "false"

    cases = phaseb_case_definitions()
    haskey(cases, case_name) || error("unsupported case $(case_name)")
    case = cases[case_name]
    models, _assets = prepare_phaseb_models(case; data_root=data_root)

    workload() = run_phaseb_case(models; threaded=threaded)
    stack, samplings = workload()
    trial_times = time_samples_ns(workload, samples)

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
        "stats" => Dict(
            "median_ns" => Int(round(median(trial_times))),
            "mean_ns" => Int(round(sum(trial_times) / length(trial_times))),
            "min_ns" => minimum(trial_times),
            "max_ns" => maximum(trial_times),
            "samples" => length(trial_times),
        ),
        "output" => summarize_output(stack, samplings),
        "policy" => "Julia CPU WFIRST Phase B HLC subset timing only; TTFx excluded; uses PreparedModel with shared cached HLC assets",
    )

    open(prefix * ".json", "w") do io
        JSON3.write(io, report)
    end
    println(JSON3.pretty(JSON3.read(JSON3.write(report))))
end

main()
