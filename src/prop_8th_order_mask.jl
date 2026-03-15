"""Apply an 8th-order occulter mask and return the generated amplitude mask."""
abstract type EighthOrderMaskShape end
struct LinearMaskShape <: EighthOrderMaskShape end
struct CircularMaskShape <: EighthOrderMaskShape end
struct EllipticalMaskShape <: EighthOrderMaskShape end

@inline _mask_shape_style(circular::Bool, elliptical::Nothing) = circular ? CircularMaskShape() : LinearMaskShape()
@inline _mask_shape_style(circular::Bool, elliptical::Real) = EllipticalMaskShape()

@inline function _mask_profile(r::T, e::T, ll::T, mm::T, plml::T) where {T<:AbstractFloat}
    return plml - prop_sinc(pi * r * e / ll)^ll + (mm / ll) * prop_sinc(pi * r * e / mm)^mm
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

function prop_8th_order_mask(
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

    x = (collect(0:(nx - 1)) .- nx / 2) .* c
    y = (collect(0:(ny - 1)) .- ny / 2) .* c

    mask = similar(wf.field, RT, ny, nx)
    shape = _mask_shape_style(circular, elliptical)

    if shape isa EllipticalMaskShape
        axis_ratio = RT(float(elliptical))
        _fill_8th_mask!(shape, mask, x, y, e, ll, mm, plml, y_axis, axis_ratio)
    else
        _fill_8th_mask!(shape, mask, x, y, e, ll, mm, plml, y_axis)
    end

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

    if wf.field isa StridedMatrix
        scratch = ensure_fft_real_scratch!(wf.workspace.fft, ny, nx)
        prop_shift_center!(scratch, mask)
        wf.field .*= scratch
    else
        wf.field .*= backend_adapt(wf.field, prop_shift_center(mask))
    end
    return mask
end
