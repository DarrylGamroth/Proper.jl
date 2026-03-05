"""Hex-segment Zernike synthesis (Mahajan & Dai JOSA A 24, 2994-3016)."""
@inline function _hex_zernike_basis(iz::Int, r::T, t::T) where {T<:AbstractFloat}
    r2 = r * r
    r3 = r * r2
    r4 = r2 * r2
    r5 = r * r4
    r6 = r3 * r3

    if iz == 1
        return one(T)
    elseif iz == 2
        return T(2) * sqrt(T(6 / 5)) * r * cos(t)
    elseif iz == 3
        return T(2) * sqrt(T(6 / 5)) * r * sin(t)
    elseif iz == 4
        return sqrt(T(5 / 43)) * (T(12) * r2 - T(5))
    elseif iz == 5
        return T(2) * sqrt(T(15 / 7)) * r2 * sin(T(2) * t)
    elseif iz == 6
        return T(2) * sqrt(T(15 / 7)) * r2 * cos(T(2) * t)
    elseif iz == 7
        return T(4) * sqrt(T(42 / 3685)) * (T(25) * r3 - T(14) * r) * sin(t)
    elseif iz == 8
        return T(4) * sqrt(T(42 / 3685)) * (T(25) * r3 - T(14) * r) * cos(t)
    elseif iz == 9
        return (T(4) / T(3)) * sqrt(T(10)) * r3 * sin(T(3) * t)
    elseif iz == 10
        return T(4) * sqrt(T(70 / 103)) * r3 * cos(T(3) * t)
    elseif iz == 11
        return (T(3) / sqrt(T(1072205))) * (T(6020) * r4 - T(5140) * r2 + T(737))
    elseif iz == 12
        return (T(30) / sqrt(T(492583))) * (T(392) * r4 - T(249) * r2) * cos(T(2) * t)
    elseif iz == 13
        return (T(30) / sqrt(T(492583))) * (T(392) * r4 - T(249) * r2) * sin(T(2) * t)
    elseif iz == 14
        return (T(10) / T(3)) * sqrt(T(7 / 99258181)) * (
            T(10) * ((T(297) - T(598) * r2) * r2 * cos(T(2) * t)) + T(5413) * r4 * cos(T(4) * t)
        )
    elseif iz == 15
        return (T(10) / T(3)) * sqrt(T(7 / 99258181)) * (
            -T(10) * ((T(297) - T(598) * r2) * r2 * sin(T(2) * t)) + T(5413) * r4 * sin(T(4) * t)
        )
    elseif iz == 16
        return T(2) * sqrt(T(6 / 1089382547)) * (T(70369) * r - T(322280) * r3 + T(309540) * r5) * cos(t)
    elseif iz == 17
        return T(2) * sqrt(T(6 / 1089382547)) * (T(70369) * r - T(322280) * r3 + T(309540) * r5) * sin(t)
    elseif iz == 18
        return T(4) * sqrt(T(385 / 295894589)) * (T(4365) * r5 - T(3322) * r3) * cos(T(3) * t)
    elseif iz == 19
        return T(4) * sqrt(T(5 / 97)) * (T(35) * r5 - T(22) * r3) * sin(T(3) * t)
    elseif iz == 20
        return ((T(-2.17600248) * r + T(13.23551876) * r3 - T(16.15533716) * r5) * cos(t) + T(5.95928883) * r5 * cos(T(5) * t))
    elseif iz == 21
        return ((T(2.17600248) * r - T(13.23551876) * r3 + T(16.15533716) * r5) * sin(t) + T(5.95928883) * r5 * sin(T(5) * t))
    elseif iz == 22
        return T(70.01749250) * r6 - T(93.07966445) * r4 + T(33.14780774) * r2 - T(2.47059083)
    else
        return zero(T)
    end
end

"""Return a summed hex-Zernike aberration map (meters)."""
function prop_hex_zernikes(
    zindex,
    zcoeff,
    n::Integer,
    dx::Real,
    hexrad::Real,
    xhex::Real=0.0,
    yhex::Real=0.0;
    rotation::Real=0.0,
    kwargs...,
)
    ang_raw = haskey(kwargs, :ROTATION) ? kwargs[:ROTATION] : rotation

    zidx = Int.(collect(zindex))
    zc_raw = collect(zcoeff)
    length(zidx) == length(zc_raw) || throw(ArgumentError("zindex and zcoeff lengths must match"))

    T = float(promote_type(eltype(zc_raw), typeof(dx), typeof(hexrad)))
    dx_t = T(dx)
    hexrad_t = T(hexrad)
    xhex_t = T(xhex)
    yhex_t = T(yhex)
    ang = T(ang_raw)

    npx = Int(n)
    hex_xc_pix = round(Int, xhex_t / dx_t)
    hex_yc_pix = round(Int, yhex_t / dx_t)
    hex_rad_pix = trunc(Int, hexrad_t / dx_t)

    rpix = hex_rad_pix + 2
    dpix = 2 * rpix + 1

    x1 = hex_xc_pix - rpix
    y1 = hex_yc_pix - rpix

    zc = convert(Vector{T}, zc_raw)
    zer = zeros(T, dpix, dpix)
    deg_to_rad = T(pi / 180)

    @inbounds for j in 1:dpix
        xx = (x1 + (j - 1)) * dx_t - xhex_t
        for i in 1:dpix
            yy = (y1 + (i - 1)) * dx_t - yhex_t
            r = hypot(xx, yy) / hexrad_t
            t = atan(yy, xx) - ang * deg_to_rad
            acc = zero(T)
            for k in eachindex(zidx)
                acc += zc[k] * _hex_zernike_basis(zidx[k], r, t)
            end
            zer[i, j] = acc
        end
    end

    abmap = zeros(T, npx, npx)

    xleft = x1 + (npx ÷ 2 + 1)
    xright = xleft + dpix - 1
    ybottom = y1 + (npx ÷ 2 + 1)
    ytop = ybottom + dpix - 1

    if xleft <= npx && ybottom <= npx && xright >= 1 && ytop >= 1
        ax1 = max(1, xleft)
        ax2 = min(npx, xright)
        ay1 = max(1, ybottom)
        ay2 = min(npx, ytop)

        zx1 = 1 + (ax1 - xleft)
        zx2 = dpix - (xright - ax2)
        zy1 = 1 + (ay1 - ybottom)
        zy2 = dpix - (ytop - ay2)

        abmap[ay1:ay2, ax1:ax2] .= zer[zy1:zy2, zx1:zx2]
    end

    return abmap
end
