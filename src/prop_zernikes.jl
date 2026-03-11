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

@inline function _obscured_zernike_mode(j::Int, r::T, t::T, eps::T) where {T<:AbstractFloat}
    r2 = r * r
    r3 = r * r2
    r4 = r2 * r2
    r5 = r * r4
    r6 = r3 * r3
    eps2 = eps * eps
    eps4 = eps2 * eps2
    eps6 = eps4 * eps2
    eps8 = eps4 * eps4
    eps10 = eps8 * eps2
    eps12 = eps6 * eps6

    if j == 1
        return one(T)
    elseif j == 2
        return (2r * cos(t)) / sqrt(1 + eps2)
    elseif j == 3
        return (2r * sin(t)) / sqrt(1 + eps2)
    elseif j == 4
        return (sqrt(T(3)) * (1 + eps2 - 2r2)) / (-1 + eps2)
    elseif j == 5
        return (sqrt(T(6)) * r2 * sin(2t)) / sqrt(1 + eps2 + eps4)
    elseif j == 6
        return (sqrt(T(6)) * r2 * cos(2t)) / sqrt(1 + eps2 + eps4)
    elseif j == 7
        return (2sqrt(T(2)) * r * (2 + 2eps4 - 3r2 + eps2 * (2 - 3r2)) * sin(t)) / ((-1 + eps2) * sqrt(1 + 5eps2 + 5eps4 + eps6))
    elseif j == 8
        return (2sqrt(T(2)) * r * (2 + 2eps4 - 3r2 + eps2 * (2 - 3r2)) * cos(t)) / ((-1 + eps2) * sqrt(1 + 5eps2 + 5eps4 + eps6))
    elseif j == 9
        return (2sqrt(T(2)) * r3 * sin(3t)) / sqrt(1 + eps2 + eps4 + eps6)
    elseif j == 10
        return (2sqrt(T(2)) * r3 * cos(3t)) / sqrt(1 + eps2 + eps4 + eps6)
    elseif j == 11
        return (sqrt(T(5)) * (1 + eps4 - 6r2 + 6r4 + eps2 * (4 - 6r2))) / ((-1 + eps2)^2)
    elseif j == 12
        return (sqrt(T(10)) * r2 * (3 + 3eps6 - 4r2 + eps2 * (3 - 4r2) + eps4 * (3 - 4r2)) * cos(2t)) /
               ((-1 + eps2) * sqrt((1 + eps2 + eps4) * (1 + 4eps2 + 10eps4 + 4eps6 + eps8)))
    elseif j == 13
        return (sqrt(T(10)) * r2 * (3 + 3eps6 - 4r2 + eps2 * (3 - 4r2) + eps4 * (3 - 4r2)) * sin(2t)) /
               ((-1 + eps2) * sqrt((1 + eps2 + eps4) * (1 + 4eps2 + 10eps4 + 4eps6 + eps8)))
    elseif j == 14
        return (sqrt(T(10)) * r4 * cos(4t)) / sqrt(1 + eps2 + eps4 + eps6 + eps8)
    elseif j == 15
        return (sqrt(T(10)) * r4 * sin(4t)) / sqrt(1 + eps2 + eps4 + eps6 + eps8)
    elseif j == 16
        return (2sqrt(T(3)) * r * (3 + 3eps8 - 12r2 + 10r4 - 12eps6 * (-1 + r2) + 2eps4 * (15 - 24r2 + 5r4) + 4eps2 * (3 - 12r2 + 10r4)) * cos(t)) /
               (((-1 + eps2)^2) * sqrt((1 + 4eps2 + eps4) * (1 + 9eps2 + 9eps4 + eps6)))
    elseif j == 17
        return (2sqrt(T(3)) * r * (3 + 3eps8 - 12r2 + 10r4 - 12eps6 * (-1 + r2) + 2eps4 * (15 - 24r2 + 5r4) + 4eps2 * (3 - 12r2 + 10r4)) * sin(t)) /
               (((-1 + eps2)^2) * sqrt((1 + 4eps2 + eps4) * (1 + 9eps2 + 9eps4 + eps6)))
    elseif j == 18
        return (2sqrt(T(3)) * r3 * (4 + 4eps8 - 5r2 + eps2 * (4 - 5r2) + eps4 * (4 - 5r2) + eps6 * (4 - 5r2)) * cos(3t)) /
               ((-1 + eps2) * sqrt((1 + eps2 + eps4 + eps6) * (1 + 4eps2 + 10eps4 + 20eps6 + 10eps8 + 4eps10 + eps12)))
    elseif j == 19
        return (2sqrt(T(3)) * r3 * (4 + 4eps8 - 5r2 + eps2 * (4 - 5r2) + eps4 * (4 - 5r2) + eps6 * (4 - 5r2)) * sin(3t)) /
               ((-1 + eps2) * sqrt((1 + eps2 + eps4 + eps6) * (1 + 4eps2 + 10eps4 + 20eps6 + 10eps8 + 4eps10 + eps12)))
    elseif j == 20
        return (2sqrt(T(3)) * r5 * cos(5t)) / sqrt(1 + eps2 + eps4 + eps6 + eps8 + eps10)
    elseif j == 21
        return (2sqrt(T(3)) * r5 * sin(5t)) / sqrt(1 + eps2 + eps4 + eps6 + eps8 + eps10)
    elseif j == 22
        return (sqrt(T(7)) * (1 + eps6 - 12r2 + 30r4 - 20r6 + eps4 * (9 - 12r2) + eps2 * (9 - 36r2 + 30r4))) / ((-1 + eps2)^3)
    end
    return zero(T)
end

function _zernike_maps(wf::WaveFront, nterms::Int; radius_m::Union{Nothing,Real}=nothing)
    nterms > 0 || throw(ArgumentError("nterms must be positive"))
    ny, nx = size(wf.field)
    x = coordinate_axis(nx, wf.sampling_m)
    y = coordinate_axis(ny, wf.sampling_m)

    beam_radius = radius_m === nothing ? prop_get_beamradius(wf) : float(radius_m)
    descs = prop_noll_zernikes(nterms)

    RT = real(eltype(wf.field))
    out = similar(wf.field, RT, ny, nx, nterms)
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
    nums = zernike_num isa Number ? [Int(zernike_num)] : Int.(collect(zernike_num))
    vals = zernike_val isa Number ? [float(zernike_val)] : float.(collect(zernike_val))
    length(nums) == length(vals) || throw(ArgumentError("zernike_num and zernike_val size mismatch"))
    eps_t = float(eps)
    if eps_t != 0 && maximum(nums) > 22
        throw(ArgumentError("PROP_ZERNIKES: Maximum index for an obscured Zernike polynomial is 22."))
    end

    RT = real(eltype(wf.field))
    dmap = similar(wf.field, RT, size(wf.field)...)
    fill!(dmap, zero(RT))

    if eps_t == 0
        nmax = maximum(nums)
        all_maps = _zernike_maps(wf, nmax; radius_m=radius)
        @inbounds for (j, coeff) in zip(nums, vals)
            dmap .+= coeff .* @view(all_maps[:, :, j])
        end
    else
        ny, nx = size(wf.field)
        x = coordinate_axis(nx, wf.sampling_m)
        y = coordinate_axis(ny, wf.sampling_m)
        beam_radius = radius === nothing ? prop_get_beamradius(wf) : float(radius)
        epsr = RT(eps_t)

        @inbounds for j in 1:nx
            xj = x[j] / beam_radius
            for i in 1:ny
                yi = y[i] / beam_radius
                r = hypot(xj, yi)
                t = atan(yi, xj)
                acc = zero(RT)
                for k in eachindex(nums)
                    acc += RT(vals[k]) * _obscured_zernike_mode(nums[k], r, t, epsr)
                end
                dmap[i, j] = acc
            end
        end
    end

    if !no_apply
        if amplitude
            wf.field .*= backend_adapt(wf.field, dmap)
        else
            wf.field .*= cis.((2pi / wf.wavelength_m) .* backend_adapt(wf.field, dmap))
        end
    end

    return dmap
end
