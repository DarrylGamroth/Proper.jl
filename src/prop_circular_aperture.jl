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

function prop_circular_aperture(wf::WaveFront, radius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    opts = EllipseOptions(real(eltype(wf.field)), kwargs)
    return _apply_shifted_ellipse!(wf, radius, radius, xc, yc, opts, false)
end
