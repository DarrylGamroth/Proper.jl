struct ErrorMapOptions{T<:AbstractFloat}
    sampling::Union{Nothing,T}
    xc_map::Union{Nothing,T}
    yc_map::Union{Nothing,T}
    rotatemap::Union{Nothing,T}
    magnify::Union{Nothing,T}
    multiply::Union{Nothing,T}
    nm::Bool
    microns::Bool
    amplitude::Bool
    wavefront::Bool
    mirror_surface::Bool

    function ErrorMapOptions{T}(
        sampling,
        xc_map,
        yc_map,
        rotatemap,
        magnify,
        multiply,
        nm::Bool,
        microns::Bool,
        amplitude::Bool,
        wavefront::Bool,
        mirror_surface::Bool,
    ) where {T<:AbstractFloat}
        return new{T}(
            _nothing_or(T, sampling),
            _nothing_or(T, xc_map),
            _nothing_or(T, yc_map),
            _nothing_or(T, rotatemap),
            _nothing_or(T, magnify),
            _nothing_or(T, multiply),
            nm,
            microns,
            amplitude,
            wavefront,
            mirror_surface,
        )
    end
end

@inline _nothing_or(::Type{T}, x::Nothing) where {T<:AbstractFloat} = nothing
@inline _nothing_or(::Type{T}, x::Real) where {T<:AbstractFloat} = T(x)

@inline function ErrorMapOptions(kwargs::Base.Iterators.Pairs)
    T = Float64
    sam = kw_lookup_float(kwargs, :SAMPLING, nothing)
    xcm = kw_lookup_float(kwargs, :XC_MAP, nothing)
    ycm = kw_lookup_float(kwargs, :YC_MAP, nothing)
    rmap = kw_lookup_float(kwargs, :ROTATEMAP, nothing)
    mag = kw_lookup_float(kwargs, :MAGNIFY, nothing)
    mul = kw_lookup_float(kwargs, :MULTIPLY, nothing)

    nmv = kw_lookup_bool(kwargs, :NM, false)
    micv = kw_lookup_bool(kwargs, :MICRONS, false)
    ampv = kw_lookup_bool(kwargs, :AMPLITUDE, false)
    wfv = kw_lookup_bool(kwargs, :WAVEFRONT, false)
    mirv = kw_lookup_bool(kwargs, :MIRROR, false) || kw_lookup_bool(kwargs, :MIRROR_SURFACE, false)

    ampv && (nmv || micv) && throw(ArgumentError("ERRORMAP: Cannot specify NM or MICRONS for an amplitude map"))
    nmv && micv && throw(ArgumentError("ERRORMAP: Cannot specify both NM and MICRONS"))

    return ErrorMapOptions{T}(sam, xcm, ycm, rmap, mag, mul, nmv, micv, ampv, wfv, mirv)
end

@inline function _errormap_maptype(opts::ErrorMapOptions)::Symbol
    return opts.amplitude ? :amplitude : (opts.mirror_surface ? :mirror_surface : :wavefront)
end

function _prop_errormap!(wf::WaveFront, filename::AbstractString, xshift::Real, yshift::Real, opts::ErrorMapOptions)
    ctx = RunContext(wf)
    dmap = prop_readmap(
        wf,
        filename,
        xshift,
        yshift;
        SAMPLING=opts.sampling,
        XC_MAP=opts.xc_map,
        YC_MAP=opts.yc_map,
    )

    if opts.rotatemap !== nothing || opts.magnify !== nothing
        dmap = shift_center_for_wavefront!(wf, dmap)
        if opts.rotatemap !== nothing
            dmap = prop_rotate(dmap, opts.rotatemap, ctx; METH="cubic", MISSING=0.0)
        end
        if opts.magnify !== nothing
            dmap = prop_magnify(dmap, opts.magnify, size(dmap, 1), ctx)
        end
        dmap = shift_center_for_wavefront!(wf, dmap; inverse=true)
    end

    if opts.microns
        dmap .*= 1e-6
    elseif opts.nm
        dmap .*= 1e-9
    end
    if opts.multiply !== nothing
        dmap .*= opts.multiply
    end

    maptype = _errormap_maptype(opts)
    if maptype === :amplitude
        wf.field .*= same_backend_style(typeof(wf.field), typeof(dmap)) ? dmap : backend_adapt(wf.field, dmap)
    else
        scale = maptype === :mirror_surface ? -4pi / wf.wavelength_m : 2pi / wf.wavelength_m
        map_backend = same_backend_style(typeof(wf.field), typeof(dmap)) ? dmap : backend_adapt(wf.field, dmap)
        wf.field .*= cis.(scale .* map_backend)
    end
    return wf
end

"""Apply map from FITS file to wavefront as phase or amplitude."""
function prop_errormap(
    wf::WaveFront,
    filename::AbstractString,
    xshift::Real=0.0,
    yshift::Real=0.0;
    kwargs...,
)
    opts = ErrorMapOptions(kwargs)
    return _prop_errormap!(wf, filename, xshift, yshift, opts)
end
