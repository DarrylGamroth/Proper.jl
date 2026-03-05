@inline function _facti(n::Int)
    n < 0 && return 0.0
    v = 1.0
    @inbounds for k in 2:n
        v *= k
    end
    return v
end

@inline function _radial_zernike(n::Int, m::Int, r)
    acc = zero(r)
    top = (n - m) ÷ 2
    @inbounds for s in 0:top
        num = (-1)^s * _facti(n - s)
        den = _facti(s) * _facti((n + m) ÷ 2 - s) * _facti((n - m) ÷ 2 - s)
        acc += (num / den) * r^(n - 2s)
    end
    return acc
end

@inline function _zernike_mode(desc, r::Real, t::Real)
    n = desc.n
    m = desc.m
    trig = desc.trig

    if n == 0 && m == 0
        return 1.0
    end

    radial = _radial_zernike(n, m, r)
    norm = m == 0 ? sqrt(n + 1.0) : sqrt(2.0 * (n + 1.0))

    if trig === :none
        return norm * radial
    elseif trig === :cos
        return norm * radial * cos(m * t)
    else
        return norm * radial * sin(m * t)
    end
end

function _zernike_maps(wf::WaveFront, nterms::Int; radius_m::Union{Nothing,Real}=nothing)
    nterms > 0 || throw(ArgumentError("nterms must be positive"))
    ny, nx = size(wf.field)
    x = coordinate_axis(nx, wf.sampling_m)
    y = coordinate_axis(ny, wf.sampling_m)

    beam_radius = radius_m === nothing ? prop_get_beamradius(wf) : float(radius_m)
    descs = prop_noll_zernikes(nterms)

    out = Array{Float64}(undef, ny, nx, nterms)
    @inbounds for j in 1:nx
        xj = x[j]
        for i in 1:ny
            yi = y[i]
            r = hypot(xj, yi) / beam_radius
            t = atan(yi, xj)
            for k in 1:nterms
                out[i, j, k] = _zernike_mode(descs[k], r, t)
            end
        end
    end
    return out
end

"""Return first `nterms` Noll-ordered Zernike mode maps."""
function prop_zernikes(wf::WaveFront, nterms::Integer)
    return _zernike_maps(wf, Int(nterms))
end

"""Apply Zernike aberrations, returning generated map (`no_apply=true` skips application)."""
function prop_zernikes(
    wf::WaveFront,
    zernike_num,
    zernike_val,
    eps::Real=0.0;
    amplitude::Bool=false,
    no_apply::Bool=false,
    radius::Union{Nothing,Real}=nothing,
)
    eps == 0 || throw(ArgumentError("Obscured Zernikes (eps != 0) are not implemented yet"))

    nums = zernike_num isa Number ? [Int(zernike_num)] : Int.(collect(zernike_num))
    vals = zernike_val isa Number ? [float(zernike_val)] : float.(collect(zernike_val))
    length(nums) == length(vals) || throw(ArgumentError("zernike_num and zernike_val size mismatch"))

    nmax = maximum(nums)
    all_maps = _zernike_maps(wf, nmax; radius_m=radius)

    dmap = zeros(Float64, size(wf.field))
    @inbounds for (j, coeff) in zip(nums, vals)
        dmap .+= coeff .* @view(all_maps[:, :, j])
    end

    if !no_apply
        if amplitude
            wf.field .*= dmap
        else
            wf.field .*= cis.((2pi / wf.wavelength_m) .* dmap)
        end
    end

    return dmap
end
