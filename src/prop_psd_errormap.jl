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

function _prop_psd_errormap!(wf::WaveFront, amp::Real, b::Real, c::Real, opts::PSDErrorMapOptions)
    n = size(wf.field, 1)
    dx = wf.sampling_m
    fpath = opts.file !== nothing && isfile(opts.file) ? opts.file : nothing
    maptype = "wavefront"

    dmap = if fpath === nothing
        dk = 1 / (n * dx)
        inclination = opts.inclination
        rotation = opts.rotation

        xk = repeat(collect(0:(n - 1))' .- (n ÷ 2), n, 1)
        yk = transpose(xk)

        # Python 3.3.4 behavior mutates xk before yk rotation and then reuses it.
        xk = xk .* cos(-rotation) .- yk .* sin(-rotation)
        yk = xk .* sin(-rotation) .+ yk .* cos(-rotation)
        yk .*= cos(-inclination)
        kpsd = sqrt.(xk .^ 2 .+ yk .^ 2) .* dk

        tpf = opts.tpf
        psd2d = if tpf
            float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ float(c))
        else
            float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ 2) .^ ((float(c) + 1) / 2)
        end
        psd2d[n ÷ 2 + 1, n ÷ 2 + 1] = 0.0
        rms_psd = sqrt(sum(psd2d)) * dk

        maxfreq = opts.max_frequency
        if maxfreq !== nothing
            kpsd[kpsd .> maxfreq] .= 0.0
            psd2d .*= kpsd
        end

        phase = 2pi .* rand(opts.rng, n, n) .- pi
        dmap_fft = fft(prop_shift_center(sqrt.(psd2d) ./ dk .* cis.(phase))) ./ length(psd2d)
        dmap_local = real(dmap_fft) ./ (n^2 * dx^2)
        rms_map = std(dmap_local)

        if !opts.rms && opts.amplitude_target === nothing
            dmap_local .*= rms_psd / (rms_map + eps())
        else
            dmap_local .*= float(amp) / (rms_map + eps())
        end
        dmap_local
    else
        prop_readmap(wf, fpath)
    end

    if !opts.no_apply
        if opts.amplitude_target !== nothing
            dmap .+= opts.amplitude_target - maximum(dmap)
            wf.field .*= dmap
            maptype = "amplitude"
        elseif opts.mirror
            wf.field .*= cis.((4pi / wf.wavelength_m) .* dmap)
            maptype = "mirror surface"
        else
            wf.field .*= cis.((2pi / wf.wavelength_m) .* dmap)
            maptype = "wavefront"
        end
    end

    dmap = prop_shift_center(dmap)

    if opts.file !== nothing && fpath === nothing
        fname = opts.file
        header = Dict{String,Any}(
            "MAPTYPE" => (maptype, " error map type"),
            "X_UNIT" => ("meters", " X-Y units"),
            "PIXSIZE" => (dx, " spacing in meters"),
            "PSD_AMP" => (amp, maptype == "amplitude" ? " PSD low frequency RMS amplitude (amp^2m^4)" : " PSD low frequency RMS amplitude (m^4)"),
            "PSD_B" => (b, " PSD correlation length (cycles/m)"),
            "PSD_C" => (c, " PSD high frequency power law"),
            "XC_PIX" => (n ÷ 2, " Center X pixel coordinate"),
            "YC_PIX" => (n ÷ 2, " Center Y pixel coordinate"),
        )
        if maptype != "amplitude"
            header["Z_UNIT"] = ("meters", " Error units")
        end
        if opts.tpf
            header["PSDTYPE"] = ("TPF", "")
        end
        if opts.max_frequency !== nothing
            header["MAXFREQ"] = (opts.max_frequency, " Maximum spatial frequency in cycles/meter")
        end
        prop_fits_write(fname, dmap; HEADER=header)
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
