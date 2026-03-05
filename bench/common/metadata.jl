module BenchMetadata

using Dates

export benchmark_metadata

function benchmark_metadata(; run_tag::String, backend::Symbol=:cpu, baseline::String="python334_patched")
    return Dict(
        "timestamp_utc" => string(now(UTC)),
        "run_tag" => run_tag,
        "baseline" => baseline,
        "backend" => String(backend),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
    )
end

end
