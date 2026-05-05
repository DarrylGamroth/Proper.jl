struct CircleOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
end

struct CircleGeometry{T<:AbstractFloat}
    xoffset_pix::T
    yoffset_pix::T
    rad_pix::T
    threshold_hi2::T
    threshold_lo2::T
    limit2::T
end

@inline function CircleOptions(::Type{T}, kwargs::Base.Iterators.Pairs) where {T<:AbstractFloat}
    return CircleOptions{T}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
    )
end

@inline CircleOptions(kwargs::Base.Iterators.Pairs) = CircleOptions(Float64, kwargs)

@inline function _apply_shifted_mask!(
    field::AbstractMatrix{<:Complex},
    mask::AbstractMatrix{<:Real};
    invert::Bool=false,
)
    ny, nx = size(field)
    if ka_mask_enabled(typeof(field), ny, nx)
        return ka_apply_shifted_mask!(field, backend_adapt(field, mask); invert=invert)
    end

    return _apply_shifted_mask_loop!(field, mask; invert=invert)
end

@inline function _apply_shifted_mask_loop!(
    field::AbstractMatrix{<:Complex},
    mask::AbstractMatrix{<:Real};
    invert::Bool=false,
)
    ny, nx = size(field)
    sy = ny ÷ 2
    sx = nx ÷ 2

    @inbounds for j in 1:nx
        js = mod1(j + sx, nx)
        for i in 1:ny
            is = mod1(i + sy, ny)
            m = mask[is, js]
            field[i, j] *= invert ? (one(m) - m) : m
        end
    end

    return field
end

@inline function circle_geometry(
    ::Type{T},
    wf::WaveFront,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::CircleOptions,
) where {T<:AbstractFloat}
    dx = T(prop_get_sampling(wf))
    beamrad_pix = T(prop_get_beamradius(wf)) / dx
    rad_pix = opts.norm ? T(radius) * beamrad_pix : T(radius) / dx
    iszero(rad_pix) && throw(ArgumentError("radius must be non-zero"))
    xoffset_pix = opts.norm ? T(xc) * beamrad_pix : T(xc) / dx
    yoffset_pix = opts.norm ? T(yc) * beamrad_pix : T(yc) / dx
    dr = inv(rad_pix)
    hi = rad_pix * (one(T) + dr)
    lo = rad_pix * (one(T) - dr)
    limit = rad_pix * T(1 + 1e-10)
    return CircleGeometry(
        xoffset_pix,
        yoffset_pix,
        rad_pix,
        hi * hi,
        lo > zero(T) ? (lo * lo) : -one(T),
        limit * limit,
    )
end

@inline function circle_mask_value(x0::T, y0::T, geom::CircleGeometry{T}, dark::Bool, nsub::Int) where {T<:AbstractFloat}
    r2 = x0 * x0 + y0 * y0
    pixval = if r2 > geom.threshold_hi2
        zero(T)
    elseif r2 <= geom.threshold_lo2
        one(T)
    else
        cnt = 0
        @inbounds for oy_i in 1:nsub
            ys = y0 + _ka_subsample_offset(oy_i, nsub, one(T))
            for ox_i in 1:nsub
                xs = x0 + _ka_subsample_offset(ox_i, nsub, one(T))
                cnt += ((xs * xs + ys * ys) <= geom.limit2)
            end
        end
        T(cnt) / T(nsub * nsub)
    end
    return dark ? (one(T) - pixval) : pixval
end

abstract type ShiftedCircleApplyExecStyle end
struct ShiftedCircleMaskExecStyle <: ShiftedCircleApplyExecStyle end
struct ShiftedCircleKAExecStyle <: ShiftedCircleApplyExecStyle end
abstract type CircleCenterExecStyle end
struct CenteredCircleStyle <: CircleCenterExecStyle end
struct ShiftedCircleStyle <: CircleCenterExecStyle end

@inline shifted_circle_apply_exec_style(::GeometryKAStyle) = ShiftedCircleKAExecStyle()
@inline shifted_circle_apply_exec_style(::GeometryKernelStyle) = ShiftedCircleMaskExecStyle()
@inline shifted_circle_apply_exec_style(::Type{A}) where {A<:AbstractArray} = shifted_circle_apply_exec_style(geometry_kernel_style(A))
@inline circle_center_exec_style(::FeatureEnabled) = CenteredCircleStyle()
@inline circle_center_exec_style(::FeatureDisabled) = ShiftedCircleStyle()
@inline circle_center_exec_style(geom::CircleGeometry) = circle_center_exec_style(feature_flag(iszero(geom.xoffset_pix) && iszero(geom.yoffset_pix)))

@inline ellipse_options(opts::CircleOptions{T}) where {T<:AbstractFloat} = EllipseOptions{T}(opts.norm, opts.dark, zero(T))

@inline function _apply_shifted_circle!(
    ::ShiftedCircleMaskExecStyle,
    wf::WaveFront,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::CircleOptions,
    invert::Bool,
)
    return _apply_shifted_ellipse!(wf, radius, radius, xc, yc, ellipse_options(opts), invert)
end

@inline function _apply_shifted_circle!(
    ::ShiftedCircleKAExecStyle,
    wf::WaveFront{T},
    radius::Real,
    xc::Real,
    yc::Real,
    opts::CircleOptions,
    invert::Bool,
) where {T<:AbstractFloat}
    geom = circle_geometry(T, wf, radius, xc, yc, opts)
    return _apply_shifted_circle!(circle_center_exec_style(geom), wf, geom, opts.dark, invert)
end

@inline function _apply_shifted_circle!(
    ::CenteredCircleStyle,
    wf::WaveFront,
    geom::CircleGeometry,
    dark::Bool,
    invert::Bool,
)
    ka_apply_centered_circle!(
        wf.field,
        geom.threshold_hi2,
        geom.threshold_lo2,
        geom.limit2;
        dark=dark,
        invert=invert,
        nsub=antialias_subsampling(),
    )
    return wf
end

@inline function _apply_shifted_circle!(
    ::ShiftedCircleStyle,
    wf::WaveFront,
    geom::CircleGeometry,
    dark::Bool,
    invert::Bool,
)
    ka_apply_shifted_circle!(
        wf.field,
        geom.xoffset_pix,
        geom.yoffset_pix,
        geom.rad_pix,
        geom.threshold_hi2,
        geom.threshold_lo2,
        geom.limit2;
        dark=dark,
        invert=invert,
        nsub=antialias_subsampling(),
    )
    return wf
end

@inline function _apply_shifted_circle!(
    wf::WaveFront,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::CircleOptions,
    invert::Bool,
)
    sty = shifted_circle_apply_exec_style(typeof(wf.field))
    return _apply_shifted_circle!(sty, wf, radius, xc, yc, opts, invert)
end

"""Multiply the current wavefront by a circular clear aperture."""
function prop_circular_aperture(
    wf::WaveFront,
    radius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    NORM=nothing,
    norm=nothing,
    DARK=nothing,
    dark=nothing,
    kwargs...,
)
    T = real(eltype(wf.field))
    normv = kw_resolve_bool(NORM, norm, false)
    darkv = kw_resolve_bool(DARK, dark, false)
    if !isempty(kwargs)
        normv = normv || kw_lookup_bool(kwargs, :NORM, false)
        darkv = darkv || kw_lookup_bool(kwargs, :DARK, false)
    end
    opts = CircleOptions{T}(normv, darkv)
    return _apply_shifted_circle!(wf, radius, xc, yc, opts, false)
end
