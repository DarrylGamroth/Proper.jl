using HdrHistogram
using Test

include(joinpath(@__DIR__, "..", "..", "common", "latency_histogram.jl"))
using .BenchLatency

mutable struct SequenceClock{T<:Integer}
    values::Vector{T}
    position::Int
end

SequenceClock(values::Vector{T}) where {T<:Integer} = SequenceClock{T}(values, 0)

function (clock::SequenceClock)()
    clock.position += 1
    return clock.values[clock.position]
end

@testset "latency configuration" begin
    config = LatencyConfig()
    @test config.samples_per_run == 1_000
    @test config.repetitions == 3
    @test config.percentiles == [50.0, 90.0, 99.0, 99.9]
    @test pkgversion(HdrHistogram) == v"0.4.0"

    for version in ("1.10", "1.12")
        manifest = read(
            joinpath(@__DIR__, "..", "Manifest-v$(version).toml"),
            String,
        )
        @test occursin(
            "repo-rev = \"$(BenchLatency.HDRHISTOGRAM_SOURCE_REVISION)\"",
            manifest,
        )
    end

    @test_throws ArgumentError LatencyConfig(samples_per_run=0)
    @test_throws ArgumentError LatencyConfig(repetitions=0)
    @test_throws ArgumentError LatencyConfig(warmup_samples=-1)
    @test_throws ArgumentError LatencyConfig(cooldown_seconds=Inf)
    @test_throws ArgumentError LatencyConfig(lowest_discernible_ns=0)
    @test_throws ArgumentError LatencyConfig(
        lowest_discernible_ns=10,
        highest_trackable_ns=19,
    )
    @test_throws ArgumentError LatencyConfig(significant_figures=6)
    @test_throws ArgumentError LatencyConfig(min_tail_observations=0)
    @test_throws ArgumentError LatencyConfig(percentiles=Float64[])
    @test_throws ArgumentError LatencyConfig(percentiles=[90.0, 50.0])
    @test_throws ArgumentError LatencyConfig(percentiles=[50.0, 50.0])
    @test_throws ArgumentError LatencyConfig(percentiles=[50.0, 101.0])
end

@testset "deterministic sample recording" begin
    config = LatencyConfig(
        samples_per_run=3,
        repetitions=1,
        warmup_samples=0,
        min_tail_observations=1,
    )
    histogram = new_histogram(config)
    clock = SequenceClock(UInt64[0, 10, 20, 30, 50, 60, 90, 100])
    calls = Ref(0)
    workload = () -> (calls[] += 1)

    duration_ns, gc = record_latency_samples!(histogram, workload, 3; clock)

    @test calls[] == 3
    @test duration_ns == 100
    @test HdrHistogram.total_count(histogram) == 3
    @test HdrHistogram.count_at_value(histogram, 10) == 1
    @test HdrHistogram.count_at_value(histogram, 20) == 1
    @test HdrHistogram.count_at_value(histogram, 30) == 1
    @test gc.allocated_bytes >= 0
    @test_throws ArgumentError record_latency_samples!(histogram, workload, 0; clock)

    backwards = SequenceClock(UInt64[10, 9, 8, 11])
    @test_throws ArgumentError record_latency_samples!(new_histogram(config), workload, 1; clock=backwards)

    short_range = LatencyConfig(
        samples_per_run=1,
        repetitions=1,
        warmup_samples=0,
        highest_trackable_ns=100,
        min_tail_observations=1,
    )
    out_of_range = SequenceClock(UInt64[0, 0, 1_000_000])
    @test_throws ArgumentError record_latency_samples!(
        new_histogram(short_range),
        workload,
        1;
        clock=out_of_range,
    )
end

@testset "percentile support follows sample count" begin
    config = LatencyConfig(
        samples_per_run=1_000,
        repetitions=1,
        warmup_samples=0,
        min_tail_observations=100,
    )
    histogram = new_histogram(config)
    append!(histogram, 1:1_000)
    report = histogram_report(histogram, config)

    @test report["samples"] == 1_000
    @test report["supported_percentiles"] == [50.0, 90.0]
    @test report["unsupported_percentiles"] == [99.0, 99.9]
    @test report["percentile_ns"]["50.0"] >= 500
    @test report["percentile_ns"]["90.0"] >= 900

    high_count_histogram = new_histogram(config)
    HdrHistogram.record_value!(high_count_histogram, 1, 100_000)
    high_count_report = histogram_report(high_count_histogram, config)
    @test high_count_report["supported_percentiles"] == [50.0, 90.0, 99.0, 99.9]
end

@testset "aggregation and histogram log roundtrip" begin
    config = LatencyConfig(
        samples_per_run=2,
        repetitions=2,
        warmup_samples=0,
        min_tail_observations=1,
    )
    first_histogram = new_histogram(config)
    second_histogram = new_histogram(config)
    append!(first_histogram, [10, 20])
    append!(second_histogram, [30, 40])
    HdrHistogram.tag!(first_histogram, "run_1")
    HdrHistogram.tag!(second_histogram, "run_2")
    gc = (
        allocated_bytes=Int64(0),
        malloc_calls=Int64(0),
        realloc_calls=Int64(0),
        pool_allocations=Int64(0),
        big_allocations=Int64(0),
        free_calls=Int64(0),
        gc_time_ns=Int64(0),
        gc_pauses=Int64(0),
        full_sweeps=Int64(0),
    )
    runs = [
        LatencyRun(first_histogram, 0, 100, gc),
        LatencyRun(second_histogram, 200, 100, gc),
    ]
    aggregate = aggregate_histograms(runs, config)
    measurement = LatencyMeasurement(runs, aggregate, nothing, 1_700_000_000_000, 300)

    @test eltype(measurement.runs) <: LatencyRun
    @test HdrHistogram.total_count(measurement.aggregate) == 4

    mktempdir() do directory
        path = write_histogram_log(
            joinpath(directory, "latency.hlog"),
            measurement,
            config;
            comment="test=true",
        )
        open(path, "r") do io
            reader = HdrHistogram.HistogramLogReader(io)
            decoded_first = HdrHistogram.next_interval_histogram(reader)
            decoded_second = HdrHistogram.next_interval_histogram(reader)
            @test decoded_first !== nothing
            @test decoded_second !== nothing
            @test HdrHistogram.tag(decoded_first) == "run_1"
            @test HdrHistogram.tag(decoded_second) == "run_2"
            @test HdrHistogram.total_count(decoded_first) == 2
            @test HdrHistogram.total_count(decoded_second) == 2
            @test min(decoded_first) == 10
            @test max(decoded_first) == 20
            @test min(decoded_second) == 30
            @test max(decoded_second) == 40
            @test HdrHistogram.next_interval_histogram(reader) === nothing
        end
    end
end

@testset "observer overhead smoke" begin
    config = LatencyConfig(
        samples_per_run=4,
        repetitions=1,
        warmup_samples=0,
        min_tail_observations=1,
    )
    overhead = measure_observer_overhead(config; samples=4)
    @test overhead.samples == 4
    @test overhead.clock_pair_and_call["samples"] == 4
    @test overhead.clock_pair_record_and_call["samples"] == 4
end
