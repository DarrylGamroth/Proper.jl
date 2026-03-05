using Random

struct PSDErrorMapOptions{T<:AbstractFloat,R}
    rng::R
    file::Union{Nothing,String}
    inclination::T
    rotation::T
    tpf::Bool
    max_frequency::Union{Nothing,T}
    rms::Bool
    no_apply::Bool
    mirror::Bool
    amplitude_target::Union{Nothing,T}
end

@inline function PSDErrorMapOptions(kwargs::Base.Iterators.Pairs)
    T = Float64
    rng = kw_lookup(kwargs, :RNG, Random.default_rng())
    fv = kw_lookup_string(kwargs, :FILE, nothing)
    inc = deg2rad(kw_lookup_float(kwargs, :INCLINATION, 0.0))
    rot = deg2rad(kw_lookup_float(kwargs, :ROTATION, 0.0))
    tpfv = kw_lookup_bool(kwargs, :TPF, false)
    maxf = kw_lookup_float(kwargs, :MAX_FREQUENCY, nothing)
    rmsv = kw_lookup_bool(kwargs, :RMS, false)
    noap = kw_lookup_bool(kwargs, :NO_APPLY, false)
    mir = kw_lookup_bool(kwargs, :MIRROR, false)

    av = kw_lookup(kwargs, :AMPLITUDE, nothing)
    atarget = av === nothing ? nothing : float(av)

    return PSDErrorMapOptions{T,typeof(rng)}(rng, fv, inc, rot, tpfv, maxf, rmsv, noap, mir, atarget)
end

@inline _psd_maptype(opts::PSDErrorMapOptions) =
    opts.amplitude_target === nothing ? (opts.mirror ? :mirror_surface : :wavefront) : :amplitude

@inline function _psd_maptype_string(maptype::Symbol)::String
    maptype === :amplitude && return "amplitude"
    maptype === :mirror_surface && return "mirror surface"
    return "wavefront"
end

@inline function _phase_from_map!(wf::WaveFront, dmap::AbstractMatrix, scale::Real)
    @inbounds @simd for idx in eachindex(wf.field, dmap)
        wf.field[idx] *= cis(scale * dmap[idx])
    end
    return wf
end

@inline function _read_psd_map(wf::WaveFront, fpath::AbstractString, ::Type{T}) where {T<:AbstractFloat}
    return Matrix{T}(prop_readmap(wf, fpath))
end

function _build_psd_map(
    wf::WaveFront,
    amp::Real,
    b::Real,
    c::Real,
    opts::PSDErrorMapOptions{T},
) where {T<:AbstractFloat}
    n = size(wf.field, 1)
    center = n ÷ 2
    center_pix = center + 1

    dx = T(wf.sampling_m)
    dk = inv(T(n) * dx)
    inv_dk = inv(dk)

    ampT = T(amp)
    inv_b = inv(T(b))
    cT = T(c)
    halfpow = (cT + one(T)) / T(2)

    rotation = opts.rotation
    inclination = opts.inclination
    cr = cos(-rotation)
    sr = sin(-rotation)
    ci = cos(-inclination)

    use_tpf = opts.tpf
    maxfreq = opts.max_frequency
    use_maxfreq = maxfreq !== nothing

    two_pi = T(2pi)
    piT = T(pi)

    spectrum = Matrix{Complex{T}}(undef, n, n)
    phase_u = rand(opts.rng, n, n)
    sum_psd = zero(T)

    @inbounds for j in 1:n
        x0 = T(j - 1 - center)
        for i in 1:n
            y0 = T(i - 1 - center)

            # Preserve upstream rotation-order quirk shared by Python and MATLAB.
            x1 = x0 * cr - y0 * sr
            y1 = (x1 * sr + y0 * cr) * ci
            k = sqrt(x1 * x1 + y1 * y1) * dk

            psd = if i == center_pix && j == center_pix
                zero(T)
            elseif use_tpf
                ampT / (one(T) + (k * inv_b)^cT)
            else
                ampT / (one(T) + (k * inv_b)^2)^halfpow
            end

            if use_maxfreq
                # Keep Python executable baseline behavior: multiply by zeroed k-grid.
                k_eff = k > maxfreq ? zero(T) : k
                psd *= k_eff
            end

            sum_psd += psd
            phase = muladd(two_pi, T(phase_u[i, j]), -piT)
            spectrum[i, j] = (sqrt(psd) * inv_dk) * cis(phase)
        end
    end

    spectrum = prop_shift_center(spectrum)
    FFTW.fft!(spectrum)
    spectrum ./= length(spectrum)

    inv_rx2 = inv(T(n * n) * dx * dx)
    dmap = Matrix{T}(undef, n, n)
    @inbounds @simd for idx in eachindex(dmap, spectrum)
        dmap[idx] = real(spectrum[idx]) * inv_rx2
    end

    rms_map = T(std(dmap))
    rms_psd = sqrt(sum_psd) * dk
    denom = rms_map + eps(T)
    scale = (!opts.rms && opts.amplitude_target === nothing) ? (rms_psd / denom) : (ampT / denom)
    dmap .*= scale

    return dmap
end

function _prop_psd_errormap!(wf::WaveFront, amp::Real, b::Real, c::Real, opts::PSDErrorMapOptions{T}) where {T<:AbstractFloat}
    n = size(wf.field, 1)
    dx = wf.sampling_m
    fpath = opts.file !== nothing && isfile(opts.file) ? opts.file : nothing

    dmap::Matrix{T} = if fpath === nothing
        _build_psd_map(wf, amp, b, c, opts)
    else
        _read_psd_map(wf, fpath, T)
    end

    maptype = _psd_maptype(opts)
    if opts.amplitude_target !== nothing
        dmap .+= opts.amplitude_target - maximum(dmap)
        !opts.no_apply && (wf.field .*= dmap)
    elseif !opts.no_apply
        scale = maptype === :mirror_surface ? 4pi / wf.wavelength_m : 2pi / wf.wavelength_m
        _phase_from_map!(wf, dmap, scale)
    end

    dmap = prop_shift_center(dmap)

    if opts.file !== nothing && fpath === nothing
        maptype_str = _psd_maptype_string(maptype)
        header = Dict{String,Any}(
            "MAPTYPE" => (maptype_str, " error map type"),
            "X_UNIT" => ("meters", " X-Y units"),
            "PIXSIZE" => (dx, " spacing in meters"),
            "PSD_AMP" => (amp, maptype === :amplitude ? " PSD low frequency RMS amplitude (amp^2m^4)" : " PSD low frequency RMS amplitude (m^4)"),
            "PSD_B" => (b, " PSD correlation length (cycles/m)"),
            "PSD_C" => (c, " PSD high frequency power law"),
            "XC_PIX" => (n ÷ 2, " Center X pixel coordinate"),
            "YC_PIX" => (n ÷ 2, " Center Y pixel coordinate"),
        )
        if maptype !== :amplitude
            header["Z_UNIT"] = ("meters", " Error units")
        end
        if opts.tpf
            header["PSDTYPE"] = ("TPF", "")
        end
        if opts.max_frequency !== nothing
            header["MAXFREQ"] = (opts.max_frequency, " Maximum spatial frequency in cycles/meter")
        end
        prop_fits_write(opts.file, dmap; HEADER=header)
    end

    return dmap
end

"""Create and optionally apply PSD-based surface/wavefront/amplitude error map."""
function prop_psd_errormap(
    wf::WaveFront,
    amp::Real,
    b::Real,
    c::Real;
    kwargs...,
)
    opts = PSDErrorMapOptions(kwargs)
    return _prop_psd_errormap!(wf, amp, b, c, opts)
end
