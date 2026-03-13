function cuda_device_label()
    try
        return CUDA.name(CUDA.device())
    catch
        return string(CUDA.device())
    end
end

function cuda_report_meta(run_tag::String; device::Union{Nothing,AbstractString}=nothing)
    meta = BenchMetadata.benchmark_metadata(run_tag=run_tag, backend=:cuda)
    payload = Dict{String,Any}("available" => true)
    if device !== nothing
        payload["device"] = device
    end
    return merge(meta, payload)
end

function skipped_cuda_report(run_tag::String, reason::AbstractString)
    return Dict(
        "meta" => merge(BenchMetadata.benchmark_metadata(run_tag=run_tag, backend=:cuda), Dict("available" => false)),
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

@inline function trial_stats(t)
    est = median(t)
    return Dict(
        "median_ns" => est.time,
        "median_allocs" => est.allocs,
        "median_bytes" => est.memory,
        "samples" => length(t.times),
    )
end

@inline cuda_sync() = CUDA.synchronize()
@inline cuda_zeros(::Type{T}, dims...) where {T} = CUDA.zeros(T, dims...)
@inline cuda_rand(::Type{T}, dims...) where {T} = CUDA.rand(T, dims...)

function cuda_wavefront_begin(diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=0.5)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = float(wavelength_m)
    d = float(diam)
    ndiam = n * float(beam_diam_fraction)
    sampling = d / ndiam
    field = CUDA.fill(complex(one(λ), zero(λ)), n, n)
    return WaveFront(field, λ, sampling, zero(λ), d)
end
