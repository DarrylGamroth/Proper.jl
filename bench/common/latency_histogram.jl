module BenchLatency

using HdrHistogram

export LatencyConfig,
       LatencyMeasurement,
       LatencyRun,
       aggregate_histograms,
       histogram_report,
       measure_latency,
       measure_observer_overhead,
       new_histogram,
       record_latency_samples!,
       write_histogram_log

const DEFAULT_PERCENTILES = Float64[50.0, 90.0, 99.0, 99.9]
const HDRHISTOGRAM_SOURCE_REVISION = "cd42987e0739c06d29d20b4341f088b85e7a620c"

struct LatencyConfig
    samples_per_run::Int
    repetitions::Int
    warmup_samples::Int
    cooldown_seconds::Float64
    lowest_discernible_ns::Int64
    highest_trackable_ns::Int64
    significant_figures::Int
    min_tail_observations::Int
    percentiles::Vector{Float64}
end

function LatencyConfig(;
    samples_per_run::Integer=1_000,
    repetitions::Integer=3,
    warmup_samples::Integer=3,
    cooldown_seconds::Real=0.0,
    lowest_discernible_ns::Integer=1,
    highest_trackable_ns::Integer=60_000_000_000,
    significant_figures::Integer=3,
    min_tail_observations::Integer=100,
    percentiles=DEFAULT_PERCENTILES,
)
    samples = Int(samples_per_run)
    repeats = Int(repetitions)
    warmups = Int(warmup_samples)
    cooldown = Float64(cooldown_seconds)
    lowest = Int64(lowest_discernible_ns)
    highest = Int64(highest_trackable_ns)
    figures = Int(significant_figures)
    tail_count = Int(min_tail_observations)
    requested_percentiles = Float64[percentile for percentile in percentiles]

    samples > 0 || throw(ArgumentError("samples_per_run must be positive"))
    repeats > 0 || throw(ArgumentError("repetitions must be positive"))
    warmups >= 0 || throw(ArgumentError("warmup_samples must be nonnegative"))
    isfinite(cooldown) && cooldown >= 0 ||
        throw(ArgumentError("cooldown_seconds must be finite and nonnegative"))
    lowest >= 1 || throw(ArgumentError("lowest_discernible_ns must be positive"))
    highest >= 2 * lowest ||
        throw(ArgumentError("highest_trackable_ns must be at least twice lowest_discernible_ns"))
    0 <= figures <= 5 || throw(ArgumentError("significant_figures must be between 0 and 5"))
    tail_count > 0 || throw(ArgumentError("min_tail_observations must be positive"))
    isempty(requested_percentiles) && throw(ArgumentError("percentiles must not be empty"))
    all(percentile -> isfinite(percentile) && 0.0 <= percentile <= 100.0,
        requested_percentiles) || throw(ArgumentError("percentiles must be finite and between 0 and 100"))
    issorted(requested_percentiles) || throw(ArgumentError("percentiles must be sorted ascending"))
    allunique(requested_percentiles) || throw(ArgumentError("percentiles must be unique"))

    return LatencyConfig(
        samples,
        repeats,
        warmups,
        cooldown,
        lowest,
        highest,
        figures,
        tail_count,
        requested_percentiles,
    )
end

struct LatencyRun{H,G}
    histogram::H
    start_offset_ns::Int64
    duration_ns::Int64
    gc::G
end

struct LatencyMeasurement{R,H,O}
    runs::Vector{R}
    aggregate::H
    observer_overhead::O
    start_time_epoch_ms::Int64
    total_duration_ns::Int64
end

@inline new_histogram(config::LatencyConfig) = HdrHistogram.Histogram(
    Int64,
    config.lowest_discernible_ns,
    config.highest_trackable_ns,
    config.significant_figures,
)

@inline function _clock_delta_ns(stop, start)
    stop >= start || throw(ArgumentError("benchmark clock moved backwards"))
    delta = stop - start
    delta <= typemax(Int64) || throw(OverflowError("benchmark clock delta exceeds Int64 nanoseconds"))
    return Int64(delta)
end

function _gc_report(diff::Base.GC_Diff)
    return (
        allocated_bytes=Int64(diff.allocd),
        malloc_calls=Int64(diff.malloc),
        realloc_calls=Int64(diff.realloc),
        pool_allocations=Int64(diff.poolalloc),
        big_allocations=Int64(diff.bigalloc),
        free_calls=Int64(diff.freecall),
        gc_time_ns=Int64(diff.total_time),
        gc_pauses=Int64(diff.pause),
        full_sweeps=Int64(diff.full_sweep),
    )
end

"""
    record_latency_samples!(histogram, workload, samples; clock=time_ns)

Record closed-loop service latency. Each new operation starts only after the
previous operation and its histogram update complete. The timed boundary is
the call to `workload`; histogram recording and report construction are
outside it. The returned duration includes the measurement loop's observation
overhead and is therefore suitable for the achieved completion rate.
"""
function record_latency_samples!(
    histogram::HdrHistogram.AbstractHistogram,
    workload::F,
    samples::Integer;
    clock::C=time_ns,
) where {F,C}
    sample_count = Int(samples)
    sample_count > 0 || throw(ArgumentError("samples must be positive"))

    loop_start = clock()
    gc_start = Base.gc_num()
    for _ in 1:sample_count
        sample_start = clock()
        workload()
        sample_stop = clock()
        elapsed_ns = _clock_delta_ns(sample_stop, sample_start)
        HdrHistogram.record_value!(histogram, elapsed_ns)
    end
    gc_diff = Base.GC_Diff(Base.gc_num(), gc_start)
    loop_stop = clock()

    return _clock_delta_ns(loop_stop, loop_start), _gc_report(gc_diff)
end

@inline function _percentile_supported(
    sample_count::Int64,
    percentile::Float64,
    minimum_tail_observations::Int,
)
    if percentile <= 50.0
        return sample_count >= 2
    elseif percentile >= 100.0
        return false
    end
    decimal_percentile = rationalize(BigInt, percentile; tol=eps(percentile))
    tail_fraction = (big(100) - decimal_percentile) / 100
    tail_observations = fld(
        big(sample_count) * numerator(tail_fraction),
        denominator(tail_fraction),
    )
    return tail_observations >= minimum_tail_observations
end

function histogram_report(histogram::HdrHistogram.AbstractHistogram, config::LatencyConfig)
    sample_count = Int64(HdrHistogram.total_count(histogram))
    supported = Float64[]
    unsupported = Float64[]
    for percentile in config.percentiles
        destination = _percentile_supported(
            sample_count,
            percentile,
            config.min_tail_observations,
        ) ? supported : unsupported
        push!(destination, percentile)
    end

    values = Vector{Int64}(undef, length(supported))
    isempty(supported) || HdrHistogram.value_at_percentile(histogram, supported, values)
    percentile_values = Dict{String,Int64}()
    for (percentile, value) in zip(supported, values)
        percentile_values[string(percentile)] = value
    end

    return Dict{String,Any}(
        "unit" => "nanoseconds",
        "samples" => sample_count,
        "lowest_discernible_ns" => config.lowest_discernible_ns,
        "highest_trackable_ns" => config.highest_trackable_ns,
        "significant_figures" => config.significant_figures,
        "min_ns" => Int64(min(histogram)),
        "max_ns" => Int64(max(histogram)),
        "mean_ns" => HdrHistogram.mean(histogram),
        "percentile_ns" => percentile_values,
        "supported_percentiles" => supported,
        "unsupported_percentiles" => unsupported,
        "minimum_tail_observations" => config.min_tail_observations,
    )
end

function aggregate_histograms(runs::AbstractVector{<:LatencyRun}, config::LatencyConfig)
    aggregate = new_histogram(config)
    for run in runs
        HdrHistogram.add!(aggregate, run.histogram)
    end
    return aggregate
end

function measure_observer_overhead(config::LatencyConfig; samples::Integer=min(config.samples_per_run, 10_000))
    sample_count = Int(samples)
    sample_count > 0 || throw(ArgumentError("observer overhead samples must be positive"))

    clock_histogram = new_histogram(config)
    record_target = new_histogram(config)
    clock_and_record_histogram = new_histogram(config)

    noop = () -> nothing
    record_one = () -> HdrHistogram.record_value!(record_target, 1)
    noop()
    record_one()
    record_latency_samples!(clock_histogram, noop, sample_count)
    record_latency_samples!(clock_and_record_histogram, record_one, sample_count)

    return (
        samples=sample_count,
        clock_pair_and_call=histogram_report(clock_histogram, config),
        clock_pair_record_and_call=histogram_report(clock_and_record_histogram, config),
    )
end

function measure_latency(
    workload::F,
    config::LatencyConfig;
    clock::C=time_ns,
    start_time_epoch_ms::Union{Nothing,Integer}=nothing,
    measure_overhead::Bool=true,
) where {F,C}
    for _ in 1:config.warmup_samples
        workload()
    end

    overhead = measure_overhead ? measure_observer_overhead(config) : nothing
    GC.gc()
    epoch_ms = start_time_epoch_ms === nothing ?
               round(Int64, time() * 1_000) : Int64(start_time_epoch_ms)
    benchmark_start = clock()

    runs = [begin
        histogram = new_histogram(config)
        GC.gc()
        run_start = clock()
        duration_ns, gc = record_latency_samples!(
            histogram,
            workload,
            config.samples_per_run;
            clock,
        )
        HdrHistogram.tag!(histogram, "run_$(run_index)")
        run = LatencyRun(
            histogram,
            _clock_delta_ns(run_start, benchmark_start),
            duration_ns,
            gc,
        )
        if run_index < config.repetitions && config.cooldown_seconds > 0
            sleep(config.cooldown_seconds)
        end
        run
    end for run_index in 1:config.repetitions]

    benchmark_stop = clock()
    aggregate = aggregate_histograms(runs, config)
    return LatencyMeasurement(
        runs,
        aggregate,
        overhead,
        epoch_ms,
        _clock_delta_ns(benchmark_stop, benchmark_start),
    )
end

function write_histogram_log(
    path::AbstractString,
    measurement::LatencyMeasurement,
    config::LatencyConfig;
    comment::AbstractString="",
)
    mkpath(dirname(path))
    open(path, "w") do io
        writer = HdrHistogram.HistogramLogWriter(io)
        HdrHistogram.output_log_format_version(writer)
        HdrHistogram.output_start_time(writer, measurement.start_time_epoch_ms)
        HdrHistogram.output_base_time(writer, measurement.start_time_epoch_ms)
        HdrHistogram.set_base_time!(writer, measurement.start_time_epoch_ms)
        HdrHistogram.output_comment(
            writer,
            "unit=ns,lowest=$(config.lowest_discernible_ns),highest=$(config.highest_trackable_ns)," *
            "significant_figures=$(config.significant_figures),samples_per_run=$(config.samples_per_run)," *
            "repetitions=$(config.repetitions),load_model=closed_loop",
        )
        isempty(comment) || HdrHistogram.output_comment(writer, comment)
        HdrHistogram.output_legend(writer)
        for run in measurement.runs
            start_seconds = run.start_offset_ns / 1.0e9
            end_seconds = (run.start_offset_ns + run.duration_ns) / 1.0e9
            HdrHistogram.output_interval_histogram(
                writer,
                start_seconds,
                end_seconds,
                run.histogram,
                1_000_000.0,
            )
        end
    end
    return path
end

end
