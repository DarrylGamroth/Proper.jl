module BenchMetadata

using Dates

export benchmark_metadata

function benchmark_metadata(; run_tag::String, compat_mode::Symbol=:python334, backend::Symbol=:cpu)
    return Dict(
        "timestamp_utc" => string(now(UTC)),
        "run_tag" => run_tag,
        "compat_mode" => String(compat_mode),
        "backend" => String(backend),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
    )
end

end
