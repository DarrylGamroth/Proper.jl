"""Apply an 8th-order occulter mask and return the generated amplitude mask."""
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

    fratio = prop_get_fratio(wf)
    wavelength = prop_get_wavelength(wf)
    sampling = prop_get_sampling(wf)

    hwhm_ld = meters ? float(hwhm) / (fratio * wavelength) : float(hwhm)
    e = 1.788 / hwhm_ld

    ny, nx = size(wf.field)
    ll = 3.0
    mm = 1.0
    plml = (ll - mm) / ll
    c = sampling / (fratio * wavelength)

    x = (collect(0:(nx - 1)) .- nx / 2) .* c
    y = (collect(0:(ny - 1)) .- ny / 2) .* c

    linear = !circular && elliptical === nothing
    RT = real(eltype(wf.field))
    mask = similar(wf.field, RT, ny, nx)

    if linear
        if y_axis
            v = y
            line = @. plml - prop_sinc(pi * v * e / ll)^ll + (mm / ll) * prop_sinc(pi * v * e / mm)^mm
            @inbounds for j in 1:nx
                mask[:, j] .= line
            end
        else
            v = x
            line = @. plml - prop_sinc(pi * v * e / ll)^ll + (mm / ll) * prop_sinc(pi * v * e / mm)^mm
            @inbounds for i in 1:ny
                mask[i, :] .= line
            end
        end
    else
        axis_ratio = elliptical === nothing ? 1.0 : float(elliptical)
        @inbounds for j in 1:nx
            xj = x[j]
            for i in 1:ny
                yy = y[i] / axis_ratio
                r = hypot(xj, yy)
                mask[i, j] = plml - prop_sinc(pi * r * e / ll)^ll + (mm / ll) * prop_sinc(pi * r * e / mm)^mm
            end
        end
    end

    # Renormalize in intensity space then convert back to amplitude.
    mask .*= mask
    mmin = minimum(mask)
    mask .-= mmin
    mmax = maximum(mask)
    if mmax > 0
        mask ./= mmax
    end
    mask .*= (float(max_transmission) - float(min_transmission))
    mask .+= float(min_transmission)
    mask .= sqrt.(mask)

    wf.field .*= prop_shift_center(mask)
    return mask
end
