module BenchmarkReporting

using BenchmarkTools
using JSON3
using Statistics

export arg_value, has_flag, parse_int_list, parse_string_list
export trial_stats, untimed_stats, measure_samples, write_json_report

function arg_value(flag::String, default=nothing)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx < length(ARGS) || error("missing value for $(flag)")
    return ARGS[idx + 1]
end

has_flag(flag::String) = any(==(flag), ARGS)

function parse_int_list(raw, default)
    raw === nothing && return collect(Int, default)
    text = String(raw)
    isempty(strip(text)) && return Int[]
    return parse.(Int, split(text, ","))
end

function parse_string_list(raw, default)
    raw === nothing && return collect(String, default)
    text = String(raw)
    isempty(strip(text)) && return String[]
    return String.(split(text, ","))
end

function measure_samples(workload, samples::Integer)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    workload()
    GC.gc()
    return run(@benchmarkable $workload() evals=1 samples=samples)
end

function trial_stats(t::BenchmarkTools.Trial)
    est_med = median(t)
    est_mean = mean(t)
    est_min = minimum(t)
    est_max = maximum(t)
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
        "samples" => length(t.times),
    )
end

function untimed_stats()
    return Dict(
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
end

function write_json_report(path::AbstractString, report)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, report)
    end
    JSON3.pretty(JSON3.read(JSON3.write(report)))
    println()
    return report
end

end
