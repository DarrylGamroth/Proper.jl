using JSON3

function load_report(path::AbstractString)
    return JSON3.read(read(path, String), Dict{String,Any})
end

function main(paths::AbstractVector{<:AbstractString})
    isempty(paths) && throw(ArgumentError("at least one thread-topology case report is required"))
    cases = load_report.(paths)
    all(case["correctness"]["equivalent"] === true for case in cases) ||
        error("cannot aggregate a numerically inequivalent topology case")
    first(cases)["correctness"]["reference_role"] == "created" ||
        error("the first thread-topology case must create the common reference")
    all(case["correctness"]["reference_role"] == "compared" for case in Iterators.drop(cases, 1)) ||
        error("subsequent thread-topology cases must compare against the common reference")

    first_meta = first(cases)["meta"]
    identity_keys = ("workload", "grid_n", "precision", "cpu_model", "cpu_threads")
    for case in cases, key in identity_keys
        case["meta"][key] == first_meta[key] || error(
            "thread-topology case metadata differs for $key",
        )
    end
    if first_meta["workload"] == "batch"
        all(case["meta"]["batch_size"] == first_meta["batch_size"] for case in cases) ||
            error("thread-topology batch sizes differ")
    end

    best_index = argmin([case["stats"]["median_ns"] for case in cases])
    best_case = cases[best_index]
    meta = copy(first_meta)
    for key in ("threads", "fftw_threads", "blas_threads", "blas_config", "run_tag")
        pop!(meta, key, nothing)
    end
    meta["run_tag"] = "thread_topology_cpu_matrix"
    meta["matrix_cases"] = length(cases)
    first_median_ns = first(cases)["stats"]["median_ns"]
    workload = String(first_meta["workload"])

    report = Dict(
        "meta" => meta,
        "policy" => "fresh process per Julia/FFTW topology; common serialized output reference; every case must pass numerical equivalence before aggregation",
        "best" => Dict(
            "julia_threads" => best_case["meta"]["threads"],
            "fftw_threads" => best_case["meta"]["fftw_threads"],
            "blas_threads" => best_case["meta"]["blas_threads"],
            "median_ns" => best_case["stats"]["median_ns"],
            "speedup_vs_first_case" => first_median_ns / best_case["stats"]["median_ns"],
        ),
        "cases" => cases,
    )

    report_path = get(
        ENV,
        "PROPER_BENCH_REPORT",
        joinpath(
            @__DIR__,
            "..",
            "..",
            "reports",
            "julia_thread_topology_cpu_$workload.json",
        ),
    )
    mkpath(dirname(report_path))
    open(report_path, "w") do io
        JSON3.write(io, report)
    end
    println(report)
    return report
end

main(ARGS)
