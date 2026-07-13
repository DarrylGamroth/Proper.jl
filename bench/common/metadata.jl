module BenchMetadata

using Dates
using FFTW
using LinearAlgebra

export benchmark_metadata

function _cpu_model()
    try
        info = Sys.cpu_info()
        if !isempty(info)
            first_info = first(info)
            if hasproperty(first_info, :model)
                return String(getproperty(first_info, :model))
            end
        end
    catch
    end
    return isdefined(Sys, :CPU_NAME) ? String(getfield(Sys, :CPU_NAME)) : "unknown"
end

function benchmark_metadata(; run_tag::String, backend::Symbol=:cpu, baseline::String="python334_patched")
    return Dict(
        "timestamp_utc" => string(now(UTC)),
        "run_tag" => run_tag,
        "baseline" => baseline,
        "backend" => String(backend),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "fftw_threads" => FFTW.get_num_threads(),
        "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => sprint(show, BLAS.get_config()),
        "cpu_threads" => Sys.CPU_THREADS,
        "cpu_model" => _cpu_model(),
        "machine" => Sys.MACHINE,
        "kernel" => Sys.KERNEL,
        "word_size" => Sys.WORD_SIZE,
        "host" => get(ENV, "HOSTNAME", "unknown"),
    )
end

end
