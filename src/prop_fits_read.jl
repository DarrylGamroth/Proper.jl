using FITSIO

"""Read FITS image. If `header=true`, return `(array, headerdict)`."""
function prop_fits_read(fname::AbstractString; header::Bool=false)
    arr = nothing
    hdr = Dict{String,Any}()
    FITS(fname, "r") do f
        hdu = f[1]
        arr = read(hdu)
        if header
            h = FITSIO.read_header(hdu)
            @inbounds for i in eachindex(h.keys)
                hdr[string(h.keys[i])] = h.values[i]
            end
        end
    end
    data = Float64.(arr)
    return header ? (data, hdr) : data
end
