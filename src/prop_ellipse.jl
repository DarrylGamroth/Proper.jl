"""Creates an image containing an antialiased filled ellipse."""
function prop_ellipse(wf::WaveFront, xradius::Real, yradius::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    nsub = antialias_subsampling()
    n = prop_get_gridsize(wf)
    dx = prop_get_sampling(wf)
    beamrad_pix = prop_get_beamradius(wf) / dx
    norm = switch_set(:NORM; kwargs...)
    dark = switch_set(:DARK; kwargs...)
    rotation = haskey(kwargs, :ROTATION) ? float(kwargs[:ROTATION]) : haskey(kwargs, :rotation) ? float(kwargs[:rotation]) : 0.0

    RT = real(eltype(wf.field))
    xf = RT
    xcenter_pix = xf(n ÷ 2)
    ycenter_pix = xf(n ÷ 2)

    xrad_pix::RT = zero(RT)
    yrad_pix::RT = zero(RT)
    if norm
        xcenter_pix += xf(xc) * xf(beamrad_pix)
        ycenter_pix += xf(yc) * xf(beamrad_pix)
        xrad_pix = xf(xradius) * xf(beamrad_pix)
        yrad_pix = xf(yradius) * xf(beamrad_pix)
    else
        xcenter_pix += xf(xc / dx)
        ycenter_pix += xf(yc / dx)
        xrad_pix = xf(xradius / dx)
        yrad_pix = xf(yradius / dx)
    end

    t = xf(deg2rad(rotation))
    sint = sin(t)
    cost = cos(t)

    xp = RT[-xrad_pix, xrad_pix, xrad_pix, -xrad_pix]
    yp = RT[-yrad_pix, -yrad_pix, yrad_pix, yrad_pix]
    xbox = xp .* cost .- yp .* sint .+ xcenter_pix
    ybox = xp .* sint .+ yp .* cost .+ ycenter_pix

    minx_pix = clamp(round(Int, minimum(xbox)) - 1, 0, n - 1)
    maxx_pix = clamp(round(Int, maximum(xbox)) + 1, 0, n - 1)
    nx = maxx_pix - minx_pix + 1
    miny_pix = clamp(round(Int, minimum(ybox)) - 1, 0, n - 1)
    maxy_pix = clamp(round(Int, maximum(ybox)) + 1, 0, n - 1)
    ny = maxy_pix - miny_pix + 1

    x_local = Matrix{RT}(undef, ny, nx)
    y_local = Matrix{RT}(undef, ny, nx)
    @inbounds for j in 1:ny
        yv = xf(j - 1 + miny_pix) - ycenter_pix
        for i in 1:nx
            x_local[j, i] = xf(i - 1 + minx_pix) - xcenter_pix
            y_local[j, i] = yv
        end
    end

    xr = @. (x_local * cost - y_local * sint) / xrad_pix
    yr = @. (x_local * sint + y_local * cost) / yrad_pix
    r = @. sqrt(xr * xr + yr * yr)

    delx = inv(xrad_pix)
    dely = inv(yrad_pix)
    drx = delx * cost - dely * sint
    dry = delx * sint + dely * cost
    dr = max(abs(drx), abs(dry))

    mask = fill(xf(-1), ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        rv = r[j, i]
        if rv > (1 + dr)
            mask[j, i] = zero(RT)
        elseif rv <= (1 - dr)
            mask[j, i] = one(RT)
        end
    end

    edge_idx = findall(==(-one(RT)), mask)
    nsubpix = xf(nsub * nsub)
    subpix_x = Matrix{RT}(undef, nsub, nsub)
    subpix_y = Matrix{RT}(undef, nsub, nsub)
    @inbounds for j in 1:nsub
        for i in 1:nsub
            subpix_x[j, i] = (xf(i - 1 - (nsub ÷ 2)) / xf(nsub)) + xf(minx_pix) - xcenter_pix
            subpix_y[j, i] = (xf(j - 1 - (nsub ÷ 2)) / xf(nsub)) + xf(miny_pix) - ycenter_pix
        end
    end

    limit = xf(1 + 1e-10)
    @inbounds for ci in eachindex(edge_idx)
        I = edge_idx[ci]
        jy = I[1] - 1
        ix = I[2] - 1
        cnt = 0
        for sj in 1:nsub
            for si in 1:nsub
                xs = subpix_x[sj, si] + xf(ix)
                ys = subpix_y[sj, si] + xf(jy)
                xsv = (xs * cost - ys * sint) / xrad_pix
                ysv = (xs * sint + ys * cost) / yrad_pix
                cnt += ((xsv * xsv + ysv * ysv) <= limit)
            end
        end
        mask[I] = xf(cnt) / nsubpix
    end

    if dark
        image = ones(RT, n, n)
        @views image[(miny_pix + 1):(miny_pix + ny), (minx_pix + 1):(minx_pix + nx)] .= 1 .- mask
        return image
    end
    image = zeros(RT, n, n)
    @views image[(miny_pix + 1):(miny_pix + ny), (minx_pix + 1):(minx_pix + nx)] .= mask
    return image
end
