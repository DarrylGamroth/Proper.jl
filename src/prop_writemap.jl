"""Write map FITS with PROPER-compatible map metadata."""
function prop_writemap(dmap::AbstractMatrix, filename::AbstractString; kwargs...)
    if !(haskey(kwargs, :RADIUS_PIX) || haskey(kwargs, :radius_pix) || haskey(kwargs, :SAMPLING) || haskey(kwargs, :sampling))
        throw(ArgumentError("Either RADIUS_PIX or SAMPLING keyword is required"))
    end

    maptype = switch_set(:AMPLITUDE; kwargs...) ? "amplitude" : switch_set(:MIRROR; kwargs...) ? "mirror" : "wavefront"
    nx = size(dmap, 2)
    ny = size(dmap, 1)

    hdr = Dict{String,Tuple{Any,String}}(
        "MAPTYPE" => (maptype, "error map type"),
        "X_UNIT" => ("meters", "X & Y units"),
        "XC_PIX" => (nx ÷ 2, "Center X pixel coordinate"),
        "YC_PIX" => (ny ÷ 2, "Center Y pixel coordinate"),
    )
    if maptype != "amplitude"
        hdr["Z_UNIT"] = ("meters", "Error units")
    end
    if haskey(kwargs, :RADIUS_PIX)
        hdr["RADPIX"] = (kwargs[:RADIUS_PIX], "beam radius in pixels")
    elseif haskey(kwargs, :radius_pix)
        hdr["RADPIX"] = (kwargs[:radius_pix], "beam radius in pixels")
    elseif haskey(kwargs, :SAMPLING)
        hdr["PIXSIZE"] = (kwargs[:SAMPLING], "spacing in meters")
    else
        hdr["PIXSIZE"] = (kwargs[:sampling], "spacing in meters")
    end

    return prop_fits_write(filename, dmap; HEADER=hdr)
end
