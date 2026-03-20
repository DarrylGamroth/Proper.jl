"""
    prop_readmap(wf, filename, xshift=0, yshift=0; kwargs...)

Read a FITS map on the host, resample it to the wavefront sampling, and return
the shifted map.

# Notes
- FITS decoding happens on the host via `FITSIO.jl`.
- The decoded map is then promoted to the backend of `wf.field` before the
  resampling/apply path continues.
- The returned array preserves the backend of `wf.field` where feasible.
"""
function prop_readmap(
    wf::WaveFront,
    filename::AbstractString,
    xshift::Real=0,
    yshift::Real=0;
    kwargs...,
)
    dmap, header = prop_fits_read_with_header(filename)
    ctx = RunContext(wf)
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

    ny, nx = size(wf.field)
    Tout = float(promote_type(real(eltype(dmap)), typeof(wf.sampling_m), typeof(pixsize), typeof(xc), typeof(yc), typeof(xshift), typeof(yshift)))
    dmap_backend = same_backend_style(typeof(wf.field), typeof(dmap)) ? dmap : backend_adapt(wf.field, dmap)
    out = similar(dmap_backend, Tout, ny, nx)
    prop_resamplemap!(out, wf, dmap_backend, pixsize, xc, yc, ctx, xshift, yshift)
    shifted = shift_center_for_wavefront!(wf, out; inverse=true)
    out === shifted || copyto!(out, shifted)
    return out
end
