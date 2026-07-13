abstract type EighthOrderMaskShape end
struct LinearMaskShape <: EighthOrderMaskShape end
struct CircularMaskShape <: EighthOrderMaskShape end
struct EllipticalMaskShape <: EighthOrderMaskShape end

const MASK_EXTREMA_BLOCK_SIZE = 256

@inline _mask_extrema_temp_length(n::Integer) = 2 * cld(n, 2 * MASK_EXTREMA_BLOCK_SIZE)

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
        @inbounds mask[i, j] = _mask_profile(r, e, ll, mm, plml)
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
        @inbounds mask[i, j] = _mask_profile(hypot(x, y), e, ll, mm, plml)
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
        @inbounds mask[i, j] = _mask_profile(hypot(x, y / axis_ratio), e, ll, mm, plml)
    end
end

@kernel function _ka_normalize_apply_8th_mask_kernel!(
    mask,
    field,
    squared_min,
    squared_max,
    min_transmission,
    max_transmission,
    sy::Int,
    sx::Int,
    ny::Int,
    nx::Int,
)
    I = @index(Global, NTuple)
    i = I[1]
    j = I[2]

    if i <= ny && j <= nx
        @uniform begin
            denominator = squared_max - squared_min
            transmission_range = max_transmission - min_transmission
        end
        @inbounds raw = mask[i, j]
        normalized = raw * raw - squared_min
        if denominator > zero(denominator)
            normalized /= denominator
        elseif iszero(denominator)
            normalized = zero(normalized)
        end
        normalized = sqrt(normalized * transmission_range + min_transmission)
        @inbounds mask[i, j] = normalized

        # Write to the shifted field index so every mask and field element is
        # touched exactly once. This is equivalent to prop_shift_center(mask)
        # for both odd and even rectangular grids.
        is = i + sy
        js = j + sx
        is = ifelse(is > ny, is - ny, is)
        js = ifelse(js > nx, js - nx, js)
        @inbounds field[is, js] *= normalized
    end
end

@inline function ka_normalize_apply_8th_mask!(
    mask::AbstractMatrix{T},
    field::AbstractMatrix{Complex{T}},
    squared_min::T,
    squared_max::T,
    min_transmission::T,
    max_transmission::T,
) where {T<:AbstractFloat}
    size(mask) == size(field) || throw(ArgumentError("mask size must match wavefront"))
    ny, nx = size(mask)
    backend = AK.get_backend(mask)
    _ka_normalize_apply_8th_mask_kernel!(backend, (16, 16))(
        mask,
        field,
        squared_min,
        squared_max,
        min_transmission,
        max_transmission,
        fld(ny, 2),
        fld(nx, 2),
        ny,
        nx;
        ndrange=(ny, nx),
    )
    return mask
end

@inline _mask_squared_extrema_map(x::T) where {T<:AbstractFloat} = complex(x * x, x * x)

@inline function _mask_squared_extrema_combine(a::Complex{T}, b::Complex{T}) where {T<:AbstractFloat}
    # Encode the minimum in the real component and the maximum in the
    # imaginary component so AcceleratedKernels can reduce both in one pass.
    return complex(min(real(a), real(b)), max(imag(a), imag(b)))
end

@inline function _mask_squared_extrema(mask::AbstractMatrix{T}, wf::WaveFront) where {T<:AbstractFloat}
    length(mask) <= 1 && return complex(zero(T), zero(T))

    # Planned FFT paths publish `workspace.fft.scratch` as the live field, so
    # mask reductions must use separately owned storage.
    required_temp = _mask_extrema_temp_length(length(mask))
    reduction_scratch = ensure_mask_reduction_scratch!(wf.workspace.mask, required_temp)
    neutral = complex(typemax(T), typemin(T))
    return AK.mapreduce(
        _mask_squared_extrema_map,
        _mask_squared_extrema_combine,
        mask;
        init=neutral,
        neutral=neutral,
        temp=reduction_scratch,
        block_size=MASK_EXTREMA_BLOCK_SIZE,
    )
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

@inline function _normalize_apply_8th_order_mask!(
    backend::BackendStyle,
    wf::WaveFront,
    mask::AbstractMatrix,
    min_transmission::Real,
    max_transmission::Real,
)
    RT = real(eltype(wf.field))
    mask .*= mask
    mmin = AK.minimum(mask)
    mask .-= mmin
    mmax = AK.maximum(mask)
    if mmax > 0
        mask ./= mmax
    end
    mask .*= (RT(float(max_transmission)) - RT(float(min_transmission)))
    mask .+= RT(float(min_transmission))
    mask .= sqrt.(mask)
    _apply_8th_order_mask!(backend, wf, mask)
    return mask
end

@inline function _normalize_apply_8th_order_mask!(
    ::Union{CUDABackend,AMDGPUBackend},
    wf::WaveFront,
    mask::AbstractMatrix{T},
    min_transmission::Real,
    max_transmission::Real,
) where {T<:AbstractFloat}
    squared_extrema = _mask_squared_extrema(mask, wf)
    ka_normalize_apply_8th_mask!(
        mask,
        wf.field,
        real(squared_extrema),
        imag(squared_extrema),
        T(float(min_transmission)),
        T(float(max_transmission)),
    )
    return mask
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

    # Renormalize in intensity space, convert back to amplitude, and apply the
    # centered-to-origin shift. Accelerators fuse this work after one paired
    # extrema reduction; the CPU path retains its established operations.
    return _normalize_apply_8th_order_mask!(backend, wf, mask, min_transmission, max_transmission)
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
