"""Read a FITS map and resample to wavefront sampling."""
function prop_readmap(
    wf::WaveFront,
    filename::AbstractString,
    xshift::Real=0,
    yshift::Real=0;
    kwargs...,
)
    dmap, header = prop_fits_read_with_header(filename)
    sampling = kw_lookup_float(kwargs, :SAMPLING, nothing)
    xc_map = kw_lookup_float(kwargs, :XC_MAP, nothing)
    yc_map = kw_lookup_float(kwargs, :YC_MAP, nothing)

    pixsize = if haskey(header, "RADPIX")
        prop_get_beamradius(wf) / float(header["RADPIX"])
    elseif sampling !== nothing
        sampling
    elseif haskey(header, "PIXSIZE")
        float(header["PIXSIZE"])
    else
        throw(ArgumentError("READMAP: No pixel scale available in header or kwargs"))
    end

    xc = xc_map !== nothing ? xc_map : size(dmap, 2) ÷ 2
    yc = yc_map !== nothing ? yc_map : size(dmap, 1) ÷ 2

    out = prop_resamplemap(wf, dmap, pixsize, xc, yc, xshift, yshift)
    return prop_shift_center(out)
end
