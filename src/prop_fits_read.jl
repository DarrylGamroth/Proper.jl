using FITSIO

const FITSHeader = Dict{String,Any}

@inline function _fits_header_dict(hdu)::FITSHeader
    hdr = Dict{String,Any}()
    h = FITSIO.read_header(hdu)
    @inbounds for i in eachindex(h.keys)
        hdr[string(h.keys[i])] = h.values[i]
    end
    return hdr
end

@inline function _prop_fits_read_data(fname::AbstractString)::AbstractArray
    FITS(fname, "r") do f
        data = read(f[1])::AbstractArray
        nd = ndims(data)
        return nd <= 1 ? data : permutedims(data, Tuple(nd:-1:1))
    end
end

@inline function _prop_fits_read_data_header(fname::AbstractString)::Tuple{AbstractArray,FITSHeader}
    FITS(fname, "r") do f
        hdu = f[1]
        img_raw = read(hdu)::AbstractArray
        nd = ndims(img_raw)
        img = nd <= 1 ? img_raw : permutedims(img_raw, Tuple(nd:-1:1))
        return img, _fits_header_dict(hdu)
    end
end

"""Read FITS image and return `(array, headerdict)` as host CPU arrays/metadata."""
@inline prop_fits_read_with_header(fname::AbstractString) = _prop_fits_read_data_header(fname)

"""Read FITS image. If `header=true`, return `(array, headerdict)`. Returned arrays are host CPU arrays."""
function prop_fits_read(fname::AbstractString; header::Bool=false)
    return header ? prop_fits_read_with_header(fname) : _prop_fits_read_data(fname)
end
