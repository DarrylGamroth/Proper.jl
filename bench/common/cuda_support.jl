const CUDA_PKGID = Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA")

@inline function _require_cuda_module()
    return Base.require(CUDA_PKGID)
end

function load_cuda_backend()
    try
        cuda_mod = Base.invokelatest(_require_cuda_module)
        if !Base.invokelatest(getproperty(cuda_mod, :functional))
            return nothing, "CUDA.functional() returned false"
        end
        Base.invokelatest(getproperty(cuda_mod, :allowscalar), false)
        return cuda_mod, nothing
    catch err
        return nothing, sprint(showerror, err)
    end
end

function cuda_device_label(CUDA)
    try
        device = Base.invokelatest(getproperty(CUDA, :device))
        return Base.invokelatest(getproperty(CUDA, :name), device)
    catch
        return string(Base.invokelatest(getproperty(CUDA, :device)))
    end
end

function cuda_report_meta(run_tag::String, CUDA=nothing)
    meta = BenchMetadata.benchmark_metadata(run_tag=run_tag, backend=:cuda)
    if CUDA === nothing
        return meta
    end
    return merge(
        meta,
        Dict(
            "device" => cuda_device_label(CUDA),
            "available" => true,
        ),
    )
end

function skipped_cuda_report(run_tag::String, reason::AbstractString)
    return Dict(
        "meta" => merge(cuda_report_meta(run_tag), Dict("available" => false)),
        "policy" => "CUDA benchmark skipped",
        "reason" => reason,
    )
end

function write_benchmark_report(path::AbstractString, report)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, report)
    end
    println(report)
    return report
end

@inline function trial_stats(t::BenchmarkTools.Trial)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

function cuda_wavefront_begin(CUDA, diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=0.5)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = float(wavelength_m)
    d = float(diam)
    ndiam = n * float(beam_diam_fraction)
    sampling = d / ndiam
    field = CUDA.fill(complex(one(λ), zero(λ)), n, n)
    return WaveFront(field, λ, sampling, zero(λ), d)
end
