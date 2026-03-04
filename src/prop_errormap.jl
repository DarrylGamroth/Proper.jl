"""Apply map from FITS file to wavefront as phase or amplitude."""
function prop_errormap(wf::WaveFront, filename::AbstractString, xshift::Real=0.0, yshift::Real=0.0; kwargs...)
    dmap = prop_readmap(wf, filename, xshift, yshift; kwargs...)
    if switch_set(:AMPLITUDE; kwargs...)
        wf.field .*= dmap
    else
        scale = switch_set(:MIRROR; kwargs...) ? 4pi / wf.wavelength_m : 2pi / wf.wavelength_m
        wf.field .*= cis.(scale .* dmap)
    end
    return wf
end
