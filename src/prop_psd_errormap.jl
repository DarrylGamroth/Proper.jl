using Random

"""Generate a seeded random error map and optionally apply it."""
function prop_psd_errormap(wf::WaveFront, amp::Real, b::Real, c::Real; kwargs...)
    rng = haskey(kwargs, :rng) ? kwargs[:rng] : Random.default_rng()
    ny, nx = size(wf.field)
    # Phase-2 placeholder: white noise scaled by amp.
    dmap = randn(rng, Float64, ny, nx) .* float(amp)
    if switch_set(:NO_APPLY; kwargs...)
        return dmap
    end
    if switch_set(:AMPLITUDE; kwargs...)
        wf.field .*= dmap
    else
        wf.field .*= cis.((2pi / wf.wavelength_m) .* dmap)
    end
    return dmap
end
