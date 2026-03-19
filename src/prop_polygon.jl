struct PolygonOptions{T<:AbstractFloat}
    norm::Bool
    dark::Bool
    rotation::T
end

@inline function PolygonOptions(kwargs::Base.Iterators.Pairs)
    return PolygonOptions{Float64}(
        kw_lookup_bool(kwargs, :NORM, false),
        kw_lookup_bool(kwargs, :DARK, false),
        kw_lookup_float(kwargs, :ROTATION, 0.0),
    )
end

@inline function _regular_polygon_vertices!(
    xv::AbstractVector{T},
    yv::AbstractVector{T},
    nsides::Integer,
    rad::T,
    x0::T,
    y0::T,
    θ0::T,
) where {T<:AbstractFloat}
    Δθ = -T(2pi) / T(nsides)
    θ = θ0
    @inbounds for k in 1:nsides
        sθ, cθ = sincos(θ)
        xv[k] = x0 + rad * cθ
        yv[k] = y0 + rad * sθ
        θ += Δθ
    end
    return xv, yv
end

function _prop_polygon!(
    image::AbstractMatrix,
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::PolygonOptions,
)
    nsides >= 3 || throw(ArgumentError("nsides must be >= 3"))

    RT = real(eltype(wf.field))
    beamr = RT(prop_get_beamradius(wf))
    rad = opts.norm ? RT(radius) * beamr : RT(radius)
    x0 = opts.norm ? RT(xc) * beamr : RT(xc)
    y0 = opts.norm ? RT(yc) * beamr : RT(yc)
    θ0 = RT(deg2rad(opts.rotation))

    xv, yv = ensure_mask_vertices!(wf.workspace.mask, nsides)
    _regular_polygon_vertices!(xv, yv, nsides, rad, x0, y0, θ0)

    return _prop_irregular_polygon(image, wf, xv, yv, IrregularPolygonOptions(opts.dark, false))
end

function _prop_polygon(
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real,
    yc::Real,
    opts::PolygonOptions,
)
    RT = real(eltype(wf.field))
    image = similar(wf.field, RT, size(wf.field)...)
    return _prop_polygon!(image, wf, nsides, radius, xc, yc, opts)
end

"""
    prop_polygon!(image, wf, nsides, radius, xc=0, yc=0; kwargs...)

Write an antialiased regular polygon mask into `image`.

# Arguments
- `image`: destination mask array
- `wf`: wavefront defining the grid and sampling
- `nsides`: number of polygon sides
- `radius`: circumscribed radius
- `xc`, `yc`: polygon center

# Keywords
- `NORM`: interpret radius and center in normalized beam-radius units
- `DARK`: return the complementary dark mask
- `ROTATION`: polygon rotation angle in degrees
"""
function prop_polygon!(
    image::AbstractMatrix,
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = PolygonOptions(kwargs)
    return _prop_polygon!(image, wf, nsides, radius, xc, yc, opts)
end

"""Return antialiased regular polygon mask."""
function prop_polygon(
    wf::WaveFront,
    nsides::Integer,
    radius::Real,
    xc::Real=0.0,
    yc::Real=0.0;
    kwargs...,
)
    opts = PolygonOptions(kwargs)
    return _prop_polygon(wf, nsides, radius, xc, yc, opts)
end
