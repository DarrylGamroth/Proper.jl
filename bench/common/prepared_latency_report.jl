using FFTW
using JSON3
using LinearAlgebra
using Proper

include(joinpath(@__DIR__, "metadata.jl"))
include(joinpath(@__DIR__, "prepared_execution_workloads.jl"))
include(joinpath(@__DIR__, "latency_histogram.jl"))
using .BenchLatency
using .BenchMetadata

const LATENCY_REPORTS_DIR = normpath(joinpath(@__DIR__, "..", "reports"))
const PROPER_REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function positive_env_int(name::AbstractString, default::Integer)
    value = parse(Int, get(ENV, name, string(default)))
    value > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

function nonnegative_env_int(name::AbstractString, default::Integer)
    value = parse(Int, get(ENV, name, string(default)))
    value >= 0 || throw(ArgumentError("$name must be nonnegative"))
    return value
end

function nonnegative_env_float(name::AbstractString, default::Real)
    value = parse(Float64, get(ENV, name, string(default)))
    isfinite(value) && value >= 0 ||
        throw(ArgumentError("$name must be finite and nonnegative"))
    return value
end

function latency_precision()
    value = lowercase(get(ENV, "PROPER_LATENCY_PRECISION", "float64"))
    value in ("float64", "fp64") && return Float64
    value in ("float32", "fp32") && return Float32
    throw(ArgumentError("PROPER_LATENCY_PRECISION must be Float64/fp64 or Float32/fp32"))
end

function latency_config_from_env()
    return LatencyConfig(
        samples_per_run=positive_env_int("PROPER_LATENCY_SAMPLES", 1_000),
        repetitions=positive_env_int("PROPER_LATENCY_REPETITIONS", 3),
        warmup_samples=nonnegative_env_int("PROPER_LATENCY_WARMUP", 3),
        cooldown_seconds=nonnegative_env_float("PROPER_LATENCY_COOLDOWN_SECONDS", 0.0),
        lowest_discernible_ns=positive_env_int("PROPER_LATENCY_LOWEST_NS", 1),
        highest_trackable_ns=positive_env_int("PROPER_LATENCY_HIGHEST_NS", 60_000_000_000),
        significant_figures=nonnegative_env_int("PROPER_LATENCY_SIGNIFICANT_FIGURES", 3),
        min_tail_observations=positive_env_int("PROPER_LATENCY_MIN_TAIL_OBSERVATIONS", 100),
    )
end

function _git_output(arguments...)
    try
        command = Cmd(vcat(
            String["git", "-C", PROPER_REPOSITORY_ROOT],
            String[string(argument) for argument in arguments],
        ))
        return readchomp(command)
    catch
        return "unknown"
    end
end

function _cpu_affinity()
    Sys.islinux() || return "unavailable"
    try
        for line in eachline("/proc/self/status")
            startswith(line, "Cpus_allowed_list:") || continue
            return strip(split(line, ':'; limit=2)[2])
        end
    catch
    end
    return "unknown"
end

function latency_environment_metadata()
    return Dict{String,Any}(
        "source_revision" => _git_output("rev-parse", "HEAD"),
        "source_status" => _git_output("status", "--short"),
        "active_project" => something(Base.active_project(), "none"),
        "julia_executable" => joinpath(Sys.BINDIR, Base.julia_exename()),
        "julia_cpu_target" => get(ENV, "JULIA_CPU_TARGET", "default"),
        "cpu_affinity" => _cpu_affinity(),
        "command" => join([Base.PROGRAM_FILE; ARGS], ' '),
    )
end

function prepare_preallocated_latency_case(
    field_factory::F,
    synchronize::S,
    ::Type{T},
    grid_n::Integer,
) where {F,S,T<:AbstractFloat}
    n = Int(grid_n)
    field = field_factory(Complex{T}, n, n)
    context = RunContext(typeof(field))
    wavefront = prop_begin!(
        field,
        T(2.4),
        T(0.55e-6);
        beam_diam_fraction=T(0.5),
        context,
    )
    output = similar(field, T, n, n)
    prepared = prepare_prescription(
        preallocated_steady_state_prescription,
        T(0.55),
        n;
        context,
    )
    prepared_run = prepare_run(
        prepared;
        activate_context=false,
        wavefront,
        output,
        run_context=context,
    )
    execute! = () -> begin
        result = prop_run(prepared_run)
        synchronize()
        return result
    end
    return (; execute!, output, wavefront, context, prepared_run)
end

@inline cpu_latency_field(::Type{T}, dimensions...) where {T} = Array{T}(undef, dimensions...)

function cpu_latency_reference(::Type{T}, grid_n::Integer) where {T<:AbstractFloat}
    prepared = prepare_steady_state_model(
        T,
        grid_n,
        cpu_prepared_context(T, grid_n);
        name=:prepared_latency_cpu_oracle,
    )
    output, sampling = prop_run(prepared)
    return Array(output), sampling
end

function numerical_comparison(
    reference::AbstractArray,
    candidate::AbstractArray;
    atol::Real,
    rtol::Real,
)
    size(candidate) == size(reference) || return Dict{String,Any}(
        "equivalent" => false,
        "reason" => "shape mismatch",
        "reference_shape" => collect(size(reference)),
        "candidate_shape" => collect(size(candidate)),
        "atol" => Float64(atol),
        "rtol" => Float64(rtol),
    )

    max_abs_error = 0.0
    difference_norm_squared = 0.0
    reference_norm_squared = 0.0
    equivalent = true
    @inbounds for index in eachindex(reference, candidate)
        reference_value = reference[index]
        candidate_value = candidate[index]
        difference = abs(candidate_value - reference_value)
        max_abs_error = max(max_abs_error, Float64(difference))
        difference_norm_squared += Float64(abs2(candidate_value - reference_value))
        reference_norm_squared += Float64(abs2(reference_value))
        if !isapprox(candidate_value, reference_value; atol, rtol)
            equivalent = false
        end
    end
    denominator = iszero(reference_norm_squared) ? 1.0 : sqrt(reference_norm_squared)
    relative_l2_error = sqrt(difference_norm_squared) / denominator
    equivalent &= isfinite(max_abs_error) && isfinite(relative_l2_error)

    return Dict{String,Any}(
        "equivalent" => equivalent,
        "max_abs_error" => max_abs_error,
        "relative_l2_error" => relative_l2_error,
        "atol" => Float64(atol),
        "rtol" => Float64(rtol),
    )
end

function verify_latency_case!(
    case,
    cpu_reference,
    reference_sampling;
    atol::Real,
    rtol::Real,
)
    output, sampling = case.execute!()
    output === case.output || error("prepared run did not return its caller-owned output")
    first_output = Array(output)
    output_again, sampling_again = case.execute!()
    output_again === case.output || error("prepared run changed its caller-owned output")
    second_output = Array(output_again)

    oracle_comparison = numerical_comparison(
        cpu_reference,
        first_output;
        atol,
        rtol,
    )
    repeat_comparison = numerical_comparison(
        first_output,
        second_output;
        atol,
        rtol,
    )
    sampling_equivalent = isapprox(sampling, reference_sampling; atol, rtol) &&
                          isapprox(sampling_again, reference_sampling; atol, rtol)
    equivalent = oracle_comparison["equivalent"] &&
                 repeat_comparison["equivalent"] &&
                 sampling_equivalent

    report = Dict{String,Any}(
        "equivalent" => equivalent,
        "cpu_oracle" => oracle_comparison,
        "repeatability" => repeat_comparison,
        "sampling_equivalent" => sampling_equivalent,
        "reference_sampling_m" => reference_sampling,
        "candidate_sampling_m" => sampling,
        "output_shape" => collect(size(first_output)),
        "output_eltype" => string(eltype(first_output)),
    )
    equivalent || error("latency workload failed CPU-oracle correctness: $(JSON3.write(report))")
    return report
end

function latency_run_report(run, config::LatencyConfig, run_index::Integer)
    return Dict{String,Any}(
        "run" => Int(run_index),
        "start_offset_ns" => run.start_offset_ns,
        "duration_ns" => run.duration_ns,
        "completion_rate_hz" => config.samples_per_run * 1.0e9 / run.duration_ns,
        "histogram" => histogram_report(run.histogram, config),
        "gc" => Dict{String,Any}(string(key) => value for (key, value) in pairs(run.gc)),
    )
end

function run_prepared_latency_report(
    backend::Symbol,
    field_factory::F,
    synchronize::S;
    device::AbstractString="",
    backend_package_version::AbstractString="",
) where {F,S}
    config = latency_config_from_env()
    precision = latency_precision()
    grid_n = positive_env_int("PROPER_LATENCY_GRID_N", PREPARED_STEADY_GRID_N)
    fftw_threads = positive_env_int("PROPER_LATENCY_FFTW_THREADS", 1)
    blas_threads = positive_env_int("PROPER_LATENCY_BLAS_THREADS", 1)
    prop_fftw_threads(fftw_threads)
    BLAS.set_num_threads(blas_threads)

    default_atol = precision === Float64 ? (backend === :cpu ? 5e-13 : 5e-10) :
                   (backend === :cpu ? 5e-6 : 5e-4)
    default_rtol = precision === Float64 ? (backend === :cpu ? 5e-13 : 5e-10) :
                   (backend === :cpu ? 5e-6 : 1e-3)
    atol = nonnegative_env_float("PROPER_LATENCY_ATOL", default_atol)
    rtol = nonnegative_env_float("PROPER_LATENCY_RTOL", default_rtol)

    cpu_reference, reference_sampling = cpu_latency_reference(precision, grid_n)
    case = prepare_preallocated_latency_case(
        field_factory,
        synchronize,
        precision,
        grid_n,
    )
    correctness = verify_latency_case!(
        case,
        cpu_reference,
        reference_sampling;
        atol,
        rtol,
    )

    measurement = measure_latency(case.execute!, config)
    precision_tag = precision === Float64 ? "fp64" : "fp32"
    stem = "julia_prepared_latency_$(backend)_$(precision_tag)"
    report_path = get(
        ENV,
        "PROPER_LATENCY_REPORT",
        joinpath(LATENCY_REPORTS_DIR, "$stem.json"),
    )
    histogram_path = get(
        ENV,
        "PROPER_LATENCY_HISTOGRAM",
        joinpath(LATENCY_REPORTS_DIR, "$stem.hlog"),
    )
    write_histogram_log(
        histogram_path,
        measurement,
        config;
        comment="backend=$(backend),precision=$(precision_tag),grid_n=$(grid_n),operation=prop_run_preallocated",
    )

    run_reports = Vector{Dict{String,Any}}(undef, length(measurement.runs))
    for index in eachindex(measurement.runs)
        run_reports[index] = latency_run_report(measurement.runs[index], config, index)
    end
    total_samples = config.samples_per_run * config.repetitions
    measured_duration_ns = sum(run -> run.duration_ns, measurement.runs)
    meta = merge(
        benchmark_metadata(
            run_tag="prepared_latency_$(backend)_$(precision_tag)",
            backend=backend,
            baseline="cpu_prepared_allocating_oracle",
        ),
        latency_environment_metadata(),
        Dict{String,Any}(
            "available" => true,
            "grid_n" => grid_n,
            "precision" => string(precision),
            "execution_surface" => "PreparedRun with caller-owned wavefront and output",
            "device" => isempty(device) ? "host" : device,
            "proper_version" => string(pkgversion(Proper)),
            "fftw_version" => string(pkgversion(FFTW)),
            "hdrhistogram_version" => string(pkgversion(BenchLatency.HdrHistogram)),
            "hdrhistogram_source_revision" => BenchLatency.HDRHISTOGRAM_SOURCE_REVISION,
            "backend_package_version" => isempty(backend_package_version) ?
                "not_applicable" : backend_package_version,
        ),
    )
    report = Dict{String,Any}(
        "meta" => meta,
        "policy" => "closed-loop warmed prepared-call service latency; histogram update excluded from each timed operation; CPU result readiness or explicit GPU synchronization included; TTFx, preparation, correctness oracle, host validation transfers, histogram export, and cooldown excluded; no coordinated-omission correction; no percentile regression gate on shared hardware",
        "contract" => Dict{String,Any}(
            "operation_start" => "immediately before prop_run(prepared_run)",
            "operation_end" => backend === :cpu ?
                "prop_run returned with the caller-owned output ready" :
                "backend synchronization completed and the caller-owned output was device-ready",
            "load_model" => "closed_loop_single_outstanding",
            "concurrency" => 1,
            "samples_per_run" => config.samples_per_run,
            "repetitions" => config.repetitions,
            "warmup_samples" => config.warmup_samples,
            "cooldown_seconds" => config.cooldown_seconds,
            "coordinated_omission_correction" => false,
            "minimum_tail_observations" => config.min_tail_observations,
            "histogram_unit" => "nanoseconds",
            "raw_histogram_log" => basename(histogram_path),
        ),
        "correctness" => correctness,
        "observer_overhead" => measurement.observer_overhead,
        "runs" => run_reports,
        "aggregate" => Dict{String,Any}(
            "histogram" => histogram_report(measurement.aggregate, config),
            "samples" => total_samples,
            "measured_duration_ns" => measured_duration_ns,
            "measurement_wall_duration_ns" => measurement.total_duration_ns,
            "completion_rate_hz" => total_samples * 1.0e9 / measured_duration_ns,
        ),
    )

    mkpath(dirname(report_path))
    open(report_path, "w") do io
        JSON3.pretty(io, report)
    end
    println(JSON3.write(report))
    return report
end
