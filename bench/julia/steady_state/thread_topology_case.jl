using BenchmarkTools
using FFTW
using JSON3
using LinearAlgebra
using Proper
using Serialization

include(joinpath(@__DIR__, "..", "..", "common", "metadata.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "wavefront_state.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "core_propagation_tail.jl"))
include(joinpath(@__DIR__, "..", "..", "common", "prepared_execution_workloads.jl"))
using .BenchMetadata

function positive_env_int(name::AbstractString, default::Integer)
    value = parse(Int, get(ENV, name, string(default)))
    value > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

@inline function trial_stats(trial::BenchmarkTools.Trial)
    estimate = median(trial)
    return Dict(
        "median_ns" => estimate.time,
        "median_allocs" => estimate.allocs,
        "median_bytes" => estimate.memory,
        "samples" => length(trial.times),
    )
end

function establish_or_compare_reference(
    path::AbstractString,
    field::AbstractArray,
    grid_n::Integer;
    workload::Symbol,
    atol::Real,
    rtol::Real,
)
    candidate = Array(field)
    if !isfile(path)
        mkpath(dirname(path))
        open(path, "w") do io
            serialize(io, (; grid_n=Int(grid_n), workload, field=candidate))
        end
        return Dict(
            "reference_role" => "created",
            "equivalent" => true,
            "max_abs_error" => 0.0,
            "relative_l2_error" => 0.0,
            "atol" => float(atol),
            "rtol" => float(rtol),
        )
    end

    reference = open(deserialize, path)
    reference.grid_n == grid_n || error(
        "thread-topology reference grid $(reference.grid_n) does not match $grid_n",
    )
    reference.workload === workload || error(
        "thread-topology reference workload $(reference.workload) does not match $workload",
    )
    size(reference.field) == size(candidate) || error("thread-topology reference shape mismatch")
    eltype(reference.field) === eltype(candidate) || error("thread-topology reference eltype mismatch")

    difference = candidate .- reference.field
    max_abs_error = maximum(abs, difference)
    reference_norm = norm(reference.field)
    relative_l2_error = norm(difference) / (iszero(reference_norm) ? one(reference_norm) : reference_norm)
    equivalent = true
    @inbounds for index in eachindex(candidate, reference.field)
        if !isapprox(candidate[index], reference.field[index]; atol=atol, rtol=rtol)
            equivalent = false
            break
        end
    end
    equivalent || error(
        "thread topology changed the numerical result " *
        "(max_abs_error=$max_abs_error, relative_l2_error=$relative_l2_error)",
    )

    return Dict(
        "reference_role" => "compared",
        "equivalent" => equivalent,
        "max_abs_error" => max_abs_error,
        "relative_l2_error" => relative_l2_error,
        "atol" => float(atol),
        "rtol" => float(rtol),
    )
end

function build_thread_topology_workload(::Val{:core}, grid_n::Integer)
    wf, ctx, snapshot = cpu_core_propagation_tail_case(Float64, grid_n)
    run! = () -> (restore_and_run_core_propagation_tail!(wf, snapshot, ctx); nothing)
    result = () -> wf.field
    return run!, result, Dict{String,Any}()
end

function build_thread_topology_workload(::Val{:batch}, grid_n::Integer)
    batch_size = positive_env_int(
        "PROPER_BENCH_BATCH_SIZE",
        length(PREPARED_BATCH_WAVELENGTHS),
    )
    wavelengths = collect(range(
        first(PREPARED_BATCH_WAVELENGTHS),
        last(PREPARED_BATCH_WAVELENGTHS);
        length=batch_size,
    ))
    stack, samplings, runs = prepare_preallocated_cpu_sweep(
        wavelengths,
        grid_n,
    )
    run! = () -> (prop_run_multi!(stack, samplings, runs); nothing)
    result = () -> stack
    return run!, result, Dict{String,Any}(
        "batch_size" => batch_size,
        "storage" => "caller-owned",
    )
end

function main()
    grid_n = positive_env_int("PROPER_BENCH_GRID_N", CORE_PROPAGATION_TAIL_GRID_N)
    samples = positive_env_int("PROPER_BENCH_SAMPLES", CORE_PROPAGATION_TAIL_SAMPLES)
    fftw_threads = positive_env_int("PROPER_BENCH_FFTW_THREADS", 1)
    blas_threads = positive_env_int("PROPER_BENCH_BLAS_THREADS", 1)
    expected_julia_threads = positive_env_int("PROPER_BENCH_JULIA_THREADS", Threads.nthreads())
    expected_julia_threads == Threads.nthreads() || error(
        "expected $expected_julia_threads Julia threads, got $(Threads.nthreads())",
    )

    reference_path = get(ENV, "PROPER_BENCH_REFERENCE", "")
    isempty(reference_path) && throw(ArgumentError("PROPER_BENCH_REFERENCE is required"))
    report_path = get(ENV, "PROPER_BENCH_REPORT", "")
    isempty(report_path) && throw(ArgumentError("PROPER_BENCH_REPORT is required"))
    atol = parse(Float64, get(ENV, "PROPER_BENCH_ATOL", "5e-13"))
    rtol = parse(Float64, get(ENV, "PROPER_BENCH_RTOL", "5e-13"))
    atol >= 0 || throw(ArgumentError("PROPER_BENCH_ATOL must be nonnegative"))
    rtol >= 0 || throw(ArgumentError("PROPER_BENCH_RTOL must be nonnegative"))
    workload = Symbol(lowercase(get(ENV, "PROPER_BENCH_WORKLOAD", "core")))
    workload in (:core, :batch) || throw(ArgumentError(
        "PROPER_BENCH_WORKLOAD must be core or batch",
    ))

    prop_fftw_threads(fftw_threads)
    BLAS.set_num_threads(blas_threads)
    run_workload!, workload_result, workload_meta =
        build_thread_topology_workload(Val(workload), grid_n)

    # Plan and compile before both the correctness comparison and timing.
    run_workload!()
    run_workload!()
    correctness = establish_or_compare_reference(
        reference_path,
        workload_result(),
        grid_n;
        workload,
        atol,
        rtol,
    )

    trial = run(@benchmarkable $run_workload!() evals=1 samples=samples)

    report = Dict(
        "meta" => merge(
            benchmark_metadata(run_tag="thread_topology_cpu_case", backend=:cpu),
            Dict(
                "grid_n" => grid_n,
                "precision" => "Float64",
                "workload" => String(workload),
            ),
            workload_meta,
        ),
        "policy" => "fresh-process Julia/FFTW topology case; warmed steady-state workload; BLAS fixed explicitly; TTFx excluded; elementwise numerical equivalence required before timing",
        "correctness" => correctness,
        "stats" => trial_stats(trial),
    )

    mkpath(dirname(report_path))
    open(report_path, "w") do io
        JSON3.write(io, report)
    end
    println(report)
    return report
end

main()
