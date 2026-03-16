@inline function _dm_kw(kwargs, s::Symbol, default)
    lower = Symbol(lowercase(String(s)))
    if haskey(kwargs, s)
        return kwargs[s]
    elseif haskey(kwargs, lower)
        return kwargs[lower]
    end
    return default
end

@inline function _fftconvolve_same(a::AbstractMatrix{T}, b::AbstractMatrix{T}) where {T<:AbstractFloat}
    ay, ax = size(a)
    by, bx = size(b)
    py = ay + by - 1
    px = ax + bx - 1

    pa = zeros(T, py, px)
    pb = zeros(T, py, px)
    @views pa[1:ay, 1:ax] .= a
    @views pb[1:by, 1:bx] .= b

    full = real.(ifft(fft(pa) .* fft(pb)))
    sy = (by ÷ 2) + 1
    sx = (bx ÷ 2) + 1
    @views return full[sy:(sy + ay - 1), sx:(sx + ax - 1)]
end

"""
Apply deformable mirror phase map in meters (wavefront mode by default).

This convenience path expects `dm_map` to already be sampled on the wavefront
grid and directly applies phase to `wf.field`.
"""
function prop_dm(wf::WaveFront, dm_map::AbstractMatrix; mirror::Bool=false)
    size(dm_map) == size(wf.field) || throw(ArgumentError("dm_map size must match wavefront"))
    scale = mirror ? -4pi / wf.wavelength_m : 2pi / wf.wavelength_m
    wf.field .*= cis.(scale .* backend_adapt(wf.field, dm_map))
    return wf
end

"""
Upstream-compatible deformable mirror model.

`dm_z0` is a DM surface map in actuator space (meters) or a FITS filename.
Returns the DM surface map reprojected onto the wavefront sampling.
"""
function prop_dm(
    wf::WaveFront,
    dm_z0::Union{AbstractString,AbstractMatrix},
    dm_xc::Real,
    dm_yc::Real,
    spacing::Real=0.0;
    kwargs...,
)
    if switch_set(:ZYX; kwargs...) && switch_set(:XYZ; kwargs...)
        throw(ArgumentError("PROP_DM: Cannot specify both XYZ and ZYX rotation orders"))
    end
    xyz = !switch_set(:ZYX; kwargs...)

    if switch_set(:FLIP_LR; kwargs...) && switch_set(:FLIP_UD; kwargs...)
        throw(ArgumentError("PROP_DM: Cannot specify both FLIP_LR and FLIP_UD"))
    end
    flip_lr = switch_set(:FLIP_LR; kwargs...)
    flip_ud = switch_set(:FLIP_UD; kwargs...)

    xtilt = float(_dm_kw(kwargs, :XTILT, 0.0))
    ytilt = float(_dm_kw(kwargs, :YTILT, 0.0))
    ztilt = float(_dm_kw(kwargs, :ZTILT, 0.0))

    influence_path = String(_dm_kw(kwargs, :INFLUENCE_FUNCTION_FILE, joinpath(@__DIR__, "..", "data", "influence_dm5v2_1.fits")))
    dm_z = dm_z0 isa AbstractString ? float.(prop_fits_read(String(dm_z0))) : float.(dm_z0)

    n = size(wf.field, 1)
    dx_surf = wf.sampling_m

    inf, hdr = prop_fits_read(influence_path; header=true)
    dx_inf_native = float(hdr["P2PD_M"])
    dx_dm_inf = float(hdr["C2CD_M"])
    inf_mag = round(Int, dx_dm_inf / dx_inf_native)

    ny_inf, nx_inf = size(inf)
    isodd(nx_inf) && isodd(ny_inf) || throw(ArgumentError("PROP_DM: Influence function width/height must be odd"))
    # Upstream cubic-convolution coordinates are 0-based.
    xc_inf = nx_inf ÷ 2
    yc_inf = ny_inf ÷ 2

    if flip_lr
        dm_z = reverse(dm_z; dims=2)
        inf = reverse(inf; dims=2)
    elseif flip_ud
        dm_z = reverse(dm_z; dims=1)
        inf = reverse(inf; dims=1)
    end

    has_nact = haskey(kwargs, :N_ACT_ACROSS_PUPIL) || haskey(kwargs, :n_act_across_pupil)
    if spacing != 0 && has_nact
        throw(ArgumentError("PROP_DM: Cannot specify both spacing and N_ACT_ACROSS_PUPIL"))
    end
    if spacing == 0 && !has_nact
        throw(ArgumentError("PROP_DM: Must specify either spacing or N_ACT_ACROSS_PUPIL"))
    end

    dx_dm = has_nact ? (2 * prop_get_beamradius(wf) / float(_dm_kw(kwargs, :N_ACT_ACROSS_PUPIL, 0))) : float(spacing)
    dx_inf = dx_inf_native * dx_dm / dx_dm_inf

    dm_cmd = if switch_set(:FIT; kwargs...)
        x = dx_dm .* (-2:2)
        inf_kernel = prop_cubic_conv(transpose(inf), x ./ dx_inf .+ xc_inf, x ./ dx_inf .+ yc_inf; grid=true)
        fitted, _ = prop_fit_dm(dm_z, inf_kernel)
        fitted
    else
        dm_z
    end

    ny_dm, nx_dm = size(dm_cmd)
    margin = 9 * inf_mag
    nx_grid = nx_dm * inf_mag + 2 * margin
    ny_grid = ny_dm * inf_mag + 2 * margin
    # Keep 0-based interpolation coordinates separate from 1-based Julia indices.
    xoff_grid0 = margin + (inf_mag ÷ 2)
    yoff_grid0 = xoff_grid0
    xoff_grid = xoff_grid0 + 1
    yoff_grid = yoff_grid0 + 1

    dm_grid = zeros(eltype(dm_cmd), ny_grid, nx_grid)
    @inbounds for iy in 1:ny_dm
        y = (iy - 1) * inf_mag + yoff_grid
        for ix in 1:nx_dm
            x = (ix - 1) * inf_mag + xoff_grid
            dm_grid[y, x] = dm_cmd[iy, ix]
        end
    end

    dm_grid = _fftconvolve_same(dm_grid, inf)

    xdim = min(round(Int, sqrt(2) * nx_grid * dx_inf / dx_surf), n)
    ydim = min(round(Int, sqrt(2) * ny_grid * dx_inf / dx_surf), n)

    xax = range(-(xdim ÷ 2) * dx_surf, step=dx_surf, length=xdim)
    yax = range(-(ydim ÷ 2) * dx_surf, step=dx_surf, length=ydim)
    x = repeat(reshape(xax, 1, :), ydim, 1)
    y = repeat(reshape(yax, :, 1), 1, xdim)

    a = deg2rad(float(xtilt))
    b = deg2rad(float(ytilt))
    g = deg2rad(float(ztilt))
    m = if xyz
        [
            cos(b) * cos(g) -cos(b) * sin(g) sin(b) 0.0
            cos(a) * sin(g) + sin(a) * sin(b) * cos(g) cos(a) * cos(g) - sin(a) * sin(b) * sin(g) -sin(a) * cos(b) 0.0
            sin(a) * sin(g) - cos(a) * sin(b) * cos(g) sin(a) * cos(g) + cos(a) * sin(b) * sin(g) cos(a) * cos(b) 0.0
            0.0 0.0 0.0 1.0
        ]
    else
        [
            cos(b) * cos(g) cos(g) * sin(a) * sin(b) - cos(a) * sin(g) cos(a) * cos(g) * sin(b) + sin(a) * sin(g) 0.0
            cos(b) * sin(g) cos(a) * cos(g) + sin(a) * sin(b) * sin(g) -cos(g) * sin(a) + cos(a) * sin(b) * sin(g) 0.0
            -sin(b) cos(b) * sin(a) cos(a) * cos(b) 0.0
            0.0 0.0 0.0 1.0
        ]
    end

    edge = [
        -1.0 -1.0 0.0 0.0
        1.0 -1.0 0.0 0.0
        1.0 1.0 0.0 0.0
        -1.0 1.0 0.0 0.0
    ]
    new_xyz = edge * m

    dx_dxs = (new_xyz[1, 1] - new_xyz[2, 1]) / (edge[1, 1] - edge[2, 1])
    dx_dys = (new_xyz[2, 1] - new_xyz[3, 1]) / (edge[2, 2] - edge[3, 2])
    dy_dxs = (new_xyz[1, 2] - new_xyz[2, 2]) / (edge[1, 1] - edge[2, 1])
    dy_dys = (new_xyz[2, 2] - new_xyz[3, 2]) / (edge[2, 2] - edge[3, 2])

    denom = 1 - (dy_dxs * dx_dys) / (dx_dxs * dy_dys)
    xs = (x ./ dx_dxs .- y .* dx_dys ./ (dx_dxs * dy_dys)) ./ denom
    ys = (y ./ dy_dys .- x .* dy_dxs ./ (dx_dxs * dy_dys)) ./ denom

    xdm = (xs .+ float(dm_xc) * dx_dm) ./ dx_inf .+ xoff_grid0
    ydm = (ys .+ float(dm_yc) * dx_dm) ./ dx_inf .+ yoff_grid0

    grid = prop_cubic_conv(transpose(dm_grid), xdm, ydm; grid=false)
    dmap = zeros(eltype(grid), n, n)

    gy, gx = size(grid)
    xmin = n ÷ 2 - xdim ÷ 2 + 1
    ymin = n ÷ 2 - ydim ÷ 2 + 1
    xmax = xmin + gx - 1
    ymax = ymin + gy - 1

    1 <= xmin <= xmax <= n || throw(ArgumentError("PROP_DM: X placement out of bounds"))
    1 <= ymin <= ymax <= n || throw(ArgumentError("PROP_DM: Y placement out of bounds"))
    @views dmap[ymin:ymax, xmin:xmax] .= grid

    if !switch_set(:NO_APPLY; kwargs...)
        # Match accepted D-0036 semantics: the projected DM map must be
        # transposed before centered-map application onto the wavefront grid.
        prop_add_phase(wf, 2 .* transpose(dmap))
    end

    return dmap
end
