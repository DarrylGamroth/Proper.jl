"""Read a FITS map and resample to wavefront sampling."""
function prop_readmap(
    wf::WaveFront,
    filename::AbstractString,
    xshift::Real=0,
    yshift::Real=0;
    SAMPLING::Union{Nothing,Real}=nothing,
    sampling::Union{Nothing,Real}=nothing,
    XC_MAP::Union{Nothing,Real}=nothing,
    xc_map::Union{Nothing,Real}=nothing,
    YC_MAP::Union{Nothing,Real}=nothing,
    yc_map::Union{Nothing,Real}=nothing,
    kwargs...,
)
    dmap, header = prop_fits_read_with_header(filename)

    pixsize = if haskey(header, "RADPIX")
        prop_get_beamradius(wf) / float(header["RADPIX"])
    elseif SAMPLING !== nothing
        float(SAMPLING)
    elseif sampling !== nothing
        float(sampling)
    elseif haskey(header, "PIXSIZE")
        float(header["PIXSIZE"])
    else
        throw(ArgumentError("READMAP: No pixel scale available in header or kwargs"))
    end

    xc = XC_MAP !== nothing ? XC_MAP : (xc_map !== nothing ? xc_map : size(dmap, 2) ÷ 2)
    yc = YC_MAP !== nothing ? YC_MAP : (yc_map !== nothing ? yc_map : size(dmap, 1) ÷ 2)

    out = prop_resamplemap(wf, dmap, pixsize, xc, yc, xshift, yshift)
    return prop_shift_center(out)
end
