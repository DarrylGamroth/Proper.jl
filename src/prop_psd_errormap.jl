using Random

@inline function _kw(kwargs, s::Symbol, default)
    return haskey(kwargs, s) ? kwargs[s] : haskey(kwargs, Symbol(lowercase(String(s)))) ? kwargs[Symbol(lowercase(String(s)))] : default
end

"""Create and optionally apply PSD-based surface/wavefront/amplitude error map."""
function prop_psd_errormap(wf::WaveFront, amp::Real, b::Real, c::Real; kwargs...)
    n = size(wf.field, 1)
    dx = wf.sampling_m
    rng = _kw(kwargs, :rng, Random.default_rng())
    fpath = nothing
    if haskey(kwargs, :FILE)
        candidate = String(kwargs[:FILE])
        isfile(candidate) && (fpath = candidate)
    elseif haskey(kwargs, :file)
        candidate = String(kwargs[:file])
        isfile(candidate) && (fpath = candidate)
    end

    maptype = "wavefront"
    dmap = if fpath === nothing
        dk = 1 / (n * dx)
        inclination = deg2rad(float(_kw(kwargs, :INCLINATION, 0.0)))
        rotation = deg2rad(float(_kw(kwargs, :ROTATION, 0.0)))

        xk = repeat(collect(0:(n - 1))' .- (n ÷ 2), n, 1)
        yk = transpose(xk)

        # Python 3.3.4 behavior mutates xk before yk rotation and then reuses it.
        xk = xk .* cos(-rotation) .- yk .* sin(-rotation)
        yk = xk .* sin(-rotation) .+ yk .* cos(-rotation)
        yk .*= cos(-inclination)
        kpsd = sqrt.(xk .^ 2 .+ yk .^ 2) .* dk

        tpf = switch_set(:TPF; kwargs...)
        psd2d = if tpf
            float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ float(c))
        else
            float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ 2) .^ ((float(c) + 1) / 2)
        end
        psd2d[n ÷ 2 + 1, n ÷ 2 + 1] = 0.0
        rms_psd = sqrt(sum(psd2d)) * dk

        maxfreq = _kw(kwargs, :MAX_FREQUENCY, nothing)
        if maxfreq !== nothing
            kpsd[kpsd .> float(maxfreq)] .= 0.0
            psd2d .*= kpsd
        end

        phase = 2pi .* rand(rng, n, n) .- pi
        dmap_fft = fft(prop_shift_center(sqrt.(psd2d) ./ dk .* cis.(phase))) ./ length(psd2d)
        dmap_local = real(dmap_fft) ./ (n^2 * dx^2)
        rms_map = std(dmap_local)

        if !switch_set(:RMS; kwargs...) && !haskey(kwargs, :AMPLITUDE) && !haskey(kwargs, :amplitude)
            dmap_local .*= rms_psd / (rms_map + eps())
        else
            dmap_local .*= float(amp) / (rms_map + eps())
        end
        dmap_local
    else
        fname = haskey(kwargs, :FILE) ? String(kwargs[:FILE]) : String(kwargs[:file])
        prop_readmap(wf, fname)
    end

    no_apply = switch_set(:NO_APPLY; kwargs...)
    if !no_apply
        if haskey(kwargs, :AMPLITUDE) || haskey(kwargs, :amplitude)
            target_amp = _kw(kwargs, :AMPLITUDE, 1.0)
            dmap .+= float(target_amp) - maximum(dmap)
            wf.field .*= dmap
            maptype = "amplitude"
        elseif switch_set(:MIRROR; kwargs...)
            wf.field .*= cis.((4pi / wf.wavelength_m) .* dmap)
            maptype = "mirror surface"
        else
            wf.field .*= cis.((2pi / wf.wavelength_m) .* dmap)
            maptype = "wavefront"
        end
    end

    dmap = prop_shift_center(dmap)

    if haskey(kwargs, :FILE) && fpath === nothing
        fname = String(kwargs[:FILE])
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
        if switch_set(:TPF; kwargs...)
            header["PSDTYPE"] = ("TPF", "")
        end
        if haskey(kwargs, :MAX_FREQUENCY)
            header["MAXFREQ"] = (float(kwargs[:MAX_FREQUENCY]), " Maximum spatial frequency in cycles/meter")
        end
        prop_fits_write(fname, dmap; HEADER=header)
    end

    return dmap
end
