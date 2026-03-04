"""Read a FITS map and resample to wavefront sampling."""
function prop_readmap(wf::WaveFront, filename::AbstractString, xshift::Real=0, yshift::Real=0; kwargs...)
    dmap, header = prop_fits_read(filename; header=true)

    pixsize = if haskey(header, "RADPIX")
        prop_get_beamradius(wf) / float(header["RADPIX"])
    elseif haskey(kwargs, :SAMPLING)
        float(kwargs[:SAMPLING])
    elseif haskey(kwargs, :sampling)
        float(kwargs[:sampling])
    elseif haskey(header, "PIXSIZE")
        float(header["PIXSIZE"])
    else
        throw(ArgumentError("READMAP: No pixel scale available in header or kwargs"))
    end

    xc = haskey(kwargs, :XC_MAP) ? kwargs[:XC_MAP] : haskey(kwargs, :xc_map) ? kwargs[:xc_map] : size(dmap, 2) ÷ 2
    yc = haskey(kwargs, :YC_MAP) ? kwargs[:YC_MAP] : haskey(kwargs, :yc_map) ? kwargs[:yc_map] : size(dmap, 1) ÷ 2

    out = prop_resamplemap(wf, dmap, pixsize, xc, yc, xshift, yshift)
    return prop_shift_center(out)
end
