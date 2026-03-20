using FITSIO

@inline _host_array_io_contract(::CPUBackend) = nothing
@inline _host_array_io_contract(::BackendStyle) =
    throw(ArgumentError("FITS write currently supports host CPU arrays only"))

"""
    prop_fits_write(fname, img; kwargs...)

Write a FITS image with optional `HEADER=Dict` metadata.

# Notes
- FITS writing is currently host-only. `img` must be a CPU-backed array.
"""
function prop_fits_write(fname::AbstractString, img::AbstractArray; kwargs...)
    _host_array_io_contract(backend_style(typeof(img)))
    header = haskey(kwargs, :HEADER) ? kwargs[:HEADER] : haskey(kwargs, :header) ? kwargs[:header] : nothing
    img_external = let nd = ndims(img)
        nd <= 1 ? img : permutedims(img, Tuple(nd:-1:1))
    end
    FITS(fname, "w") do f
        write(f, img_external)
        if header !== nothing
            hdu = f[1]
            for (k, v) in header
                if v isa Tuple && length(v) >= 1
                    FITSIO.write_key(hdu, String(k), v[1], length(v) > 1 ? String(v[2]) : "")
                else
                    FITSIO.write_key(hdu, String(k), v, "")
                end
            end
        end
    end
    return nothing
end
