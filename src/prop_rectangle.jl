"""Return an antialiased filled rectangle mask."""
function prop_rectangle(wf::WaveFront, xsize::Real, ysize::Real, xc::Real=0.0, yc::Real=0.0; kwargs...)
    n = prop_get_gridsize(wf)
    dx = prop_get_sampling(wf)
    beamrad = prop_get_beamradius(wf)
    pr = beamrad / dx

    norm = switch_set(:NORM; kwargs...)
    dark = switch_set(:DARK; kwargs...)
    rotation = haskey(kwargs, :ROTATION) ? float(kwargs[:ROTATION]) : haskey(kwargs, :rotation) ? float(kwargs[:rotation]) : 0.0
    θ = deg2rad(rotation)
    cθ = cos(θ)
    sθ = sin(θ)

    xcp = n ÷ 2 + (norm ? float(xc) * pr : float(xc) / dx)
    ycp = n ÷ 2 + (norm ? float(yc) * pr : float(yc) / dx)
    xrp = (norm ? float(xsize) * pr : float(xsize) / dx) / 2
    yrp = (norm ? float(ysize) * pr : float(ysize) / dx) / 2

    RT = real(eltype(wf.field))
    image = zeros(RT, n, n)

    # Bounding box from rotated corners.
    xp = (-xrp, -xrp, xrp, xrp)
    yp = (-yrp, yrp, yrp, -yrp)
    xbox = ntuple(i -> xp[i] * cθ - yp[i] * sθ + xcp, 4)
    ybox = ntuple(i -> xp[i] * sθ + yp[i] * cθ + ycp, 4)

    minx = max(0, floor(Int, minimum(xbox) - 1))
    maxx = min(n - 1, ceil(Int, maximum(xbox) + 1))
    miny = max(0, floor(Int, minimum(ybox) - 1))
    maxy = min(n - 1, ceil(Int, maximum(ybox) + 1))

    nsub = antialias_subsampling()
    inv_sub = 1 / (nsub * nsub)

    @inbounds for ypix in miny:maxy
        y0 = (ypix - ycp)
        for xpix in minx:maxx
            x0 = (xpix - xcp)
            cnt = 0
            for ys in 1:nsub
                yo = y0 + (ys - (nsub + 1) / 2) / nsub
                for xs in 1:nsub
                    xo = x0 + (xs - (nsub + 1) / 2) / nsub
                    xr = xo * cθ - yo * sθ
                    yr = xo * sθ + yo * cθ
                    cnt += (abs(xr) <= xrp && abs(yr) <= yrp)
                end
            end
            image[ypix + 1, xpix + 1] = cnt * inv_sub
        end
    end

    if dark
        image .= 1 .- image
    end

    return image
end
