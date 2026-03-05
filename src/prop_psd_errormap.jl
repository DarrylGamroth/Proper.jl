using Random

@inline function _kw(kwargs, s::Symbol, default)
    return haskey(kwargs, s) ? kwargs[s] : haskey(kwargs, Symbol(lowercase(String(s)))) ? kwargs[Symbol(lowercase(String(s)))] : default
end

"""Create and optionally apply PSD-based surface/wavefront/amplitude error map."""
function prop_psd_errormap(wf::WaveFront, amp::Real, b::Real, c::Real; kwargs...)
    n = size(wf.field, 1)
    dx = wf.sampling_m
    dk = 1 / (n * dx)

    rng = _kw(kwargs, :rng, Random.default_rng())
    inclination = deg2rad(float(_kw(kwargs, :INCLINATION, 0.0)))
    rotation = deg2rad(float(_kw(kwargs, :ROTATION, 0.0)))

    kx = repeat(collect(0:(n - 1))' .- (n ÷ 2), n, 1)
    ky = repeat(collect(0:(n - 1)) .- (n ÷ 2), 1, n)

    xr = kx .* cos(-rotation) .- ky .* sin(-rotation)
    yr = kx .* sin(-rotation) .+ ky .* cos(-rotation)
    yr .*= cos(-inclination)

    kpsd = sqrt.(xr .^ 2 .+ yr .^ 2) .* dk

    tpf = switch_set(:TPF; kwargs...)
    psd2d = if tpf
        float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ float(c))
    else
        float(amp) ./ (1 .+ (kpsd ./ float(b)) .^ 2) .^ ((float(c) + 1) / 2)
    end

    psd2d[n ÷ 2 + 1, n ÷ 2 + 1] = 0.0

    maxfreq = _kw(kwargs, :MAX_FREQUENCY, nothing)
    if maxfreq !== nothing
        psd2d .*= (kpsd .<= float(maxfreq))
    end

    phase = 2pi .* rand(rng, n, n) .- pi
    dmap_fft = ifftshift(sqrt.(psd2d) ./ dk .* cis.(phase))
    dmap = real(ifft(dmap_fft) .* length(dmap_fft)) ./ (n^2 * dx^2)

    rms_map = std(dmap)
    if switch_set(:RMS; kwargs...) || switch_set(:AMPLITUDE; kwargs...)
        dmap .*= float(amp) / (rms_map + eps())
    else
        rms_psd = sqrt(sum(psd2d)) * dk
        dmap .*= rms_psd / (rms_map + eps())
    end

    no_apply = switch_set(:NO_APPLY; kwargs...)
    if !no_apply
        if switch_set(:AMPLITUDE; kwargs...)
            target_amp = _kw(kwargs, :AMPLITUDE, 1.0)
            dmap .+= float(target_amp) - maximum(dmap)
            wf.field .*= dmap
        elseif switch_set(:MIRROR; kwargs...)
            wf.field .*= cis.((4pi / wf.wavelength_m) .* dmap)
        else
            wf.field .*= cis.((2pi / wf.wavelength_m) .* dmap)
        end
    end

    return prop_shift_center(dmap)
end
