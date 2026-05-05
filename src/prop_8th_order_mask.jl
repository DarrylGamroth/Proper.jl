abstract type EighthOrderMaskShape end
struct LinearMaskShape <: EighthOrderMaskShape end
struct CircularMaskShape <: EighthOrderMaskShape end
struct EllipticalMaskShape <: EighthOrderMaskShape end

@inline _mask_shape_style(circular::Bool, elliptical::Nothing) = circular ? CircularMaskShape() : LinearMaskShape()
@inline _mask_shape_style(circular::Bool, elliptical::Real) = EllipticalMaskShape()

@inline function _mask_profile(r::T, e::T, ll::T, mm::T, plml::T) where {T<:AbstractFloat}
    return plml - prop_sinc(pi * r * e / ll)^ll + (mm / ll) * prop_sinc(pi * r * e / mm)^mm
end

@inline function _mask_axis(n::Integer, c::T) where {T<:AbstractFloat}
    center = fld(n, 2) + 1
    return range((1 - center) * c, step=c, length=n)
end

@inline function _mask_coord(idx::Int, n::Int, c::T) where {T<:AbstractFloat}
    return (T(idx) - T(fld(n, 2) + 1)) * c
end

@kernel function _ka_fill_8th_linear_mask_kernel!(
    mask,
    e,
    ll,
    mm,
    plml,
    c,
    y_axis::Bool,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        r = y_axis ? _mask_coord(i, ny, c) : _mask_coord(j, nx, c)
        mask[i, j] = _mask_profile(r, e, ll, mm, plml)
    end
end

@kernel function _ka_fill_8th_circular_mask_kernel!(
    mask,
    e,
    ll,
    mm,
    plml,
    c,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        x = _mask_coord(j, nx, c)
        y = _mask_coord(i, ny, c)
        mask[i, j] = _mask_profile(hypot(x, y), e, ll, mm, plml)
    end
end

@kernel function _ka_fill_8th_elliptical_mask_kernel!(
    mask,
    e,
    ll,
    mm,
    plml,
    c,
    axis_ratio,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        x = _mask_coord(j, nx, c)
        y = _mask_coord(i, ny, c)
        mask[i, j] = _mask_profile(hypot(x, y / axis_ratio), e, ll, mm, plml)
    end
end

function _fill_8th_mask!(::LinearMaskShape, mask::AbstractMatrix{T}, x::AbstractVector{T}, y::AbstractVector{T}, e::T, ll::T, mm::T, plml::T, y_axis::Bool) where {T<:AbstractFloat}
    if y_axis
        line = similar(y)
        @inbounds for i in eachindex(y)
            line[i] = _mask_profile(y[i], e, ll, mm, plml)
        end
        @inbounds for j in axes(mask, 2)
            mask[:, j] .= line
        end
    else
        line = similar(x)
        @inbounds for j in eachindex(x)
            line[j] = _mask_profile(x[j], e, ll, mm, plml)
        end
        @inbounds for i in axes(mask, 1)
            mask[i, :] .= line
        end
    end
    return mask
end

function _fill_8th_mask!(::CircularMaskShape, mask::AbstractMatrix{T}, x::AbstractVector{T}, y::AbstractVector{T}, e::T, ll::T, mm::T, plml::T, y_axis::Bool) where {T<:AbstractFloat}
    @inbounds for j in eachindex(x)
        xj = x[j]
        for i in eachindex(y)
            r = hypot(xj, y[i])
            mask[i, j] = _mask_profile(r, e, ll, mm, plml)
        end
    end
    return mask
end

function _fill_8th_mask!(::EllipticalMaskShape, mask::AbstractMatrix{T}, x::AbstractVector{T}, y::AbstractVector{T}, e::T, ll::T, mm::T, plml::T, y_axis::Bool, axis_ratio::T) where {T<:AbstractFloat}
    @inbounds for j in eachindex(x)
        xj = x[j]
        for i in eachindex(y)
            r = hypot(xj, y[i] / axis_ratio)
            mask[i, j] = _mask_profile(r, e, ll, mm, plml)
        end
    end
    return mask
end

@inline function _fill_8th_mask_for_shape!(
    ::EllipticalMaskShape,
    mask::AbstractMatrix{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    y_axis::Bool,
    elliptical::Real,
) where {T<:AbstractFloat}
    axis_ratio = T(float(elliptical))
    return _fill_8th_mask!(EllipticalMaskShape(), mask, x, y, e, ll, mm, plml, y_axis, axis_ratio)
end

@inline function _fill_8th_mask_for_shape!(
    shape::Union{LinearMaskShape,CircularMaskShape},
    mask::AbstractMatrix{T},
    x::AbstractVector{T},
    y::AbstractVector{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    y_axis::Bool,
    elliptical::Union{Nothing,Real},
) where {T<:AbstractFloat}
    _ = elliptical
    return _fill_8th_mask!(shape, mask, x, y, e, ll, mm, plml, y_axis)
end

@inline function _fill_8th_mask_for_shape!(
    ::CPUBackend,
    shape::EighthOrderMaskShape,
    mask::AbstractMatrix{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    c::T,
    y_axis::Bool,
    elliptical::Union{Nothing,Real},
) where {T<:AbstractFloat}
    ny, nx = size(mask)
    x = _mask_axis(nx, c)
    y = _mask_axis(ny, c)
    return _fill_8th_mask_for_shape!(shape, mask, x, y, e, ll, mm, plml, y_axis, elliptical)
end

@inline function _fill_8th_mask_for_shape!(
    ::Union{CUDABackend,AMDGPUBackend},
    ::LinearMaskShape,
    mask::AbstractMatrix{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    c::T,
    y_axis::Bool,
    elliptical::Union{Nothing,Real},
) where {T<:AbstractFloat}
    _ = elliptical
    ny, nx = size(mask)
    backend = AK.get_backend(mask)
    _ka_fill_8th_linear_mask_kernel!(backend, (16, 16))(mask, e, ll, mm, plml, c, y_axis, ny, nx; ndrange=(ny, nx))
    return mask
end

@inline function _fill_8th_mask_for_shape!(
    ::Union{CUDABackend,AMDGPUBackend},
    ::CircularMaskShape,
    mask::AbstractMatrix{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    c::T,
    y_axis::Bool,
    elliptical::Union{Nothing,Real},
) where {T<:AbstractFloat}
    _ = y_axis
    _ = elliptical
    ny, nx = size(mask)
    backend = AK.get_backend(mask)
    _ka_fill_8th_circular_mask_kernel!(backend, (16, 16))(mask, e, ll, mm, plml, c, ny, nx; ndrange=(ny, nx))
    return mask
end

@inline function _fill_8th_mask_for_shape!(
    ::Union{CUDABackend,AMDGPUBackend},
    ::EllipticalMaskShape,
    mask::AbstractMatrix{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    c::T,
    y_axis::Bool,
    elliptical::Real,
) where {T<:AbstractFloat}
    _ = y_axis
    ny, nx = size(mask)
    axis_ratio = T(float(elliptical))
    backend = AK.get_backend(mask)
    _ka_fill_8th_elliptical_mask_kernel!(backend, (16, 16))(mask, e, ll, mm, plml, c, axis_ratio, ny, nx; ndrange=(ny, nx))
    return mask
end

@inline function _fill_8th_mask_for_shape!(
    backend::BackendStyle,
    shape::EighthOrderMaskShape,
    mask::AbstractMatrix{T},
    e::T,
    ll::T,
    mm::T,
    plml::T,
    c::T,
    y_axis::Bool,
    elliptical::Union{Nothing,Real},
) where {T<:AbstractFloat}
    _ = backend
    ny, nx = size(mask)
    x = _mask_axis(nx, c)
    y = _mask_axis(ny, c)
    return _fill_8th_mask_for_shape!(shape, mask, x, y, e, ll, mm, plml, y_axis, elliptical)
end

@inline function _apply_8th_order_mask!(::CPUBackend, wf::WaveFront, mask::AbstractMatrix{T}) where {T<:AbstractFloat}
    ny, nx = size(mask)
    scratch = ensure_fft_real_scratch!(wf.workspace.fft, ny, nx)
    prop_shift_center!(scratch, mask)
    wf.field .*= scratch
    return wf.field
end

@inline function _apply_8th_order_mask!(::BackendStyle, wf::WaveFront, mask::AbstractMatrix)
    wf.field .*= backend_adapt(wf.field, prop_shift_center(mask))
    return wf.field
end

"""
    prop_8th_order_mask!(mask, wf, hwhm; kwargs...)

Apply an 8th-order occulter mask to `wf`, writing the generated amplitude mask
to caller-owned `mask`.

# Arguments
- `mask`: real output buffer with the same size and backend as `wf.field`
- `wf`: wavefront to modify in place
- `hwhm`: half-width at half-maximum of the mask profile

# Keywords
- `min_transmission`, `max_transmission`: transmission range after
  normalization
- `meters`: interpret `hwhm` in meters instead of `lambda/D`
- `circular`: use the circular variant of the mask
- `elliptical`: axis ratio for the elliptical variant
- `y_axis`: apply the linear mask profile along the `y` axis instead of `x`

# Returns
- `mask`, after it has been multiplied into the wavefront.
"""
function prop_8th_order_mask!(
    mask::AbstractMatrix,
    wf::WaveFront,
    hwhm::Real;
    min_transmission::Real=0.0,
    max_transmission::Real=1.0,
    meters::Bool=false,
    circular::Bool=false,
    elliptical::Union{Nothing,Real}=nothing,
    y_axis::Bool=false,
    kwargs...,
)
    meters = meters || switch_set(:METERS; kwargs...)
    circular = circular || switch_set(:CIRCULAR; kwargs...)
    y_axis = y_axis || switch_set(:Y_AXIS; kwargs...)
    if haskey(kwargs, :ELLIPTICAL)
        elliptical = float(kwargs[:ELLIPTICAL])
    elseif haskey(kwargs, :elliptical)
        elliptical = float(kwargs[:elliptical])
    end

    RT = real(eltype(wf.field))
    size(mask) == size(wf.field) || throw(ArgumentError("mask size must match wavefront"))
    eltype(mask) === RT || throw(ArgumentError("mask eltype must be $(RT)"))
    same_backend_style(typeof(mask), typeof(wf.field)) ||
        throw(ArgumentError("mask and wavefront must use the same backend"))

    fratio = RT(prop_get_fratio(wf))
    wavelength = RT(prop_get_wavelength(wf))
    sampling = RT(prop_get_sampling(wf))

    hwhm_ld = meters ? RT(float(hwhm) / (fratio * wavelength)) : RT(float(hwhm))
    e = RT(1.788) / hwhm_ld

    ny, nx = size(wf.field)
    ll = RT(3)
    mm = RT(1)
    plml = (ll - mm) / ll
    c = sampling / (fratio * wavelength)

    shape = _mask_shape_style(circular, elliptical)
    backend = backend_style(typeof(wf.field))
    _fill_8th_mask_for_shape!(backend, shape, mask, e, ll, mm, plml, c, y_axis, elliptical)

    # Renormalize in intensity space then convert back to amplitude.
    mask .*= mask
    mmin = minimum(mask)
    mask .-= mmin
    mmax = maximum(mask)
    if mmax > 0
        mask ./= mmax
    end
    mask .*= (RT(float(max_transmission)) - RT(float(min_transmission)))
    mask .+= RT(float(min_transmission))
    mask .= sqrt.(mask)

    _apply_8th_order_mask!(backend, wf, mask)
    return mask
end

"""
    prop_8th_order_mask(wf, hwhm; kwargs...)

Apply an 8th-order occulter mask to `wf` and return the generated amplitude
mask.

# Arguments
- `wf`: wavefront to modify in place
- `hwhm`: half-width at half-maximum of the mask profile

# Keywords
- `min_transmission`, `max_transmission`: transmission range after
  normalization
- `meters`: interpret `hwhm` in meters instead of `lambda/D`
- `circular`: use the circular variant of the mask
- `elliptical`: axis ratio for the elliptical variant
- `y_axis`: apply the linear mask profile along the `y` axis instead of `x`

# Returns
- The generated amplitude mask before it is multiplied into the wavefront.
"""
function prop_8th_order_mask(wf::WaveFront, hwhm::Real; kwargs...)
    RT = real(eltype(wf.field))
    mask = similar(wf.field, RT, size(wf.field)...)
    return prop_8th_order_mask!(mask, wf, hwhm; kwargs...)
end
