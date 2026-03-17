function amdgpu_device_label()
    try
        return string(AMDGPU.device())
    catch
        return "AMDGPU device"
    end
end

function amdgpu_report_meta(run_tag::String; device::Union{Nothing,AbstractString}=nothing)
    meta = BenchMetadata.benchmark_metadata(run_tag=run_tag, backend=:amdgpu)
    payload = Dict{String,Any}("available" => true)
    if device !== nothing
        payload["device"] = device
    end
    return merge(meta, payload)
end

function skipped_amdgpu_report(run_tag::String, reason::AbstractString)
    return Dict(
        "meta" => merge(BenchMetadata.benchmark_metadata(run_tag=run_tag, backend=:amdgpu), Dict("available" => false)),
        "policy" => "AMDGPU benchmark skipped",
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

@inline amdgpu_sync() = AMDGPU.synchronize()
@inline amdgpu_zeros(::Type{T}, dims...) where {T} = AMDGPU.zeros(T, dims...)
@inline amdgpu_rand(::Type{T}, dims...) where {T} = AMDGPU.rand(T, dims...)

function amdgpu_wavefront_begin(::Type{T}, diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=0.5) where {T<:AbstractFloat}
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    λ = T(wavelength_m)
    d = T(diam)
    ndiam = T(n) * T(beam_diam_fraction)
    sampling = d / ndiam
    field = AMDGPU.fill(complex(one(λ), zero(λ)), n, n)
    return WaveFront(field, λ, sampling, zero(λ), d)
end

function amdgpu_wavefront_begin(diam::Real, wavelength_m::Real, gridsize::Integer; beam_diam_fraction::Real=0.5)
    T = float(typeof(wavelength_m))
    return amdgpu_wavefront_begin(T, diam, wavelength_m, gridsize; beam_diam_fraction=beam_diam_fraction)
end
