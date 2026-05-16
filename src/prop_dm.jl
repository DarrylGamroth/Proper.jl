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

const _DM_DIRECT_CONV_MAX_OPS = 40_000_000

@inline function _dm_active_count(a::AbstractMatrix)
    nactive = 0
    @inbounds for v in a
        nactive += !iszero(v)
    end
    return nactive
end

@inline _dm_direct_convolution_cost(nactive::Integer, inf::AbstractMatrix) = nactive * length(inf)

@inline function _dm_use_direct_convolution(nactive::Integer, inf::AbstractMatrix)
    return _dm_direct_convolution_cost(nactive, inf) <= _DM_DIRECT_CONV_MAX_OPS
end

@inline _dm_float_matrix(::Type{T}, a::AbstractMatrix{T}) where {T<:AbstractFloat} = a
@inline _dm_float_matrix(::Type{T}, a::AbstractMatrix) where {T<:AbstractFloat} = T.(a)

function _dm_influence_direct!(
    dm_grid::StridedMatrix{T},
    dm_cmd::AbstractMatrix,
    inf::AbstractMatrix{T},
    inf_mag::Integer,
    xoff_grid::Integer,
    yoff_grid::Integer,
) where {T<:AbstractFloat}
    fill!(dm_grid, zero(T))
    ny_dm, nx_dm = size(dm_cmd)
    ny_inf, nx_inf = size(inf)
    cy = ny_inf ÷ 2 + 1
    cx = nx_inf ÷ 2 + 1
    patches_inbounds =
        xoff_grid - cx + 1 >= 1 &&
        yoff_grid - cy + 1 >= 1 &&
        (nx_dm - 1) * inf_mag + xoff_grid + nx_inf - cx <= size(dm_grid, 2) &&
        (ny_dm - 1) * inf_mag + yoff_grid + ny_inf - cy <= size(dm_grid, 1)

    if patches_inbounds
        @inbounds for ix in 1:nx_dm
            x0 = (ix - 1) * inf_mag + xoff_grid
            for iy in 1:ny_dm
                coeff = T(dm_cmd[iy, ix])
                iszero(coeff) && continue
                y0 = (iy - 1) * inf_mag + yoff_grid
                for kx in 1:nx_inf
                    x = x0 + kx - cx
                    for ky in 1:ny_inf
                        y = y0 + ky - cy
                        dm_grid[y, x] += coeff * inf[ky, kx]
                    end
                end
            end
        end
    else
        @inbounds for ix in 1:nx_dm
            x0 = (ix - 1) * inf_mag + xoff_grid
            for iy in 1:ny_dm
                coeff = T(dm_cmd[iy, ix])
                iszero(coeff) && continue
                y0 = (iy - 1) * inf_mag + yoff_grid
                for kx in 1:nx_inf
                    x = x0 + kx - cx
                    1 <= x <= size(dm_grid, 2) || continue
                    for ky in 1:ny_inf
                        y = y0 + ky - cy
                        1 <= y <= size(dm_grid, 1) || continue
                        dm_grid[y, x] += coeff * inf[ky, kx]
                    end
                end
            end
        end
    end
    return dm_grid
end

function _dm_influence_grid(
    dm_cmd::AbstractMatrix,
    inf::AbstractMatrix,
    inf_mag::Integer,
    nx_grid::Integer,
    ny_grid::Integer,
    xoff_grid::Integer,
    yoff_grid::Integer,
)
    T = float(promote_type(eltype(dm_cmd), eltype(inf)))
    dm_grid = zeros(T, ny_grid, nx_grid)
    inf_grid = _dm_float_matrix(T, inf)
    nactive = _dm_active_count(dm_cmd)
    if _dm_use_direct_convolution(nactive, inf)
        return _dm_influence_direct!(dm_grid, dm_cmd, inf_grid, inf_mag, xoff_grid, yoff_grid)
    end

    @inbounds for iy in axes(dm_cmd, 1)
        y = (iy - 1) * inf_mag + yoff_grid
        for ix in axes(dm_cmd, 2)
            x = (ix - 1) * inf_mag + xoff_grid
            dm_grid[y, x] = T(dm_cmd[iy, ix])
        end
    end
    return _fftconvolve_same(dm_grid, inf_grid)
end

function _dm_cubic_conv_transposed_grid!(
    out::StridedMatrix,
    source::StridedMatrix,
    xcoords::AbstractVector,
    ycoords::AbstractVector,
)
    size(out) == (length(ycoords), length(xcoords)) || throw(ArgumentError("output size mismatch for DM interpolation"))
    source_t = transpose(source)
    @inbounds for j in eachindex(xcoords)
        x = xcoords[j]
        for i in eachindex(ycoords)
            out[i, j] = libcconv(source_t, ycoords[i], x)
        end
    end
    return out
end

function _dm_cubic_conv_transposed_coordinate_grid!(
    out::StridedMatrix,
    source::StridedMatrix,
    xgrid::AbstractMatrix,
    ygrid::AbstractMatrix,
)
    size(xgrid) == size(ygrid) || throw(ArgumentError("xgrid and ygrid sizes must match"))
    size(out) == size(xgrid) || throw(ArgumentError("output size mismatch for DM interpolation"))
    source_t = transpose(source)
    @inbounds for j in axes(xgrid, 2)
        for i in axes(xgrid, 1)
            out[i, j] = libcconv(source_t, ygrid[i, j], xgrid[i, j])
        end
    end
    return out
end

function _dm_apply_projected_phase!(wf::WaveFront, dmap::AbstractMatrix)
    return prop_add_phase(wf, 2 .* transpose(dmap))
end

function _dm_apply_projected_phase!(
    wf::WaveFront{T,<:StridedMatrix{Complex{T}}},
    dmap::StridedMatrix,
) where {T<:AbstractFloat}
    field = wf.field
    size(dmap) == size(field) || throw(ArgumentError("phase size must match wavefront"))
    nrow, ncol = size(field)
    sy = _center_shift_amount(nrow, true)
    sx = _center_shift_amount(ncol, true)
    scale = T(4pi) / wf.wavelength_m
    @inbounds for j in axes(field, 2)
        src_j = mod1(j - sx, ncol)
        for i in axes(field, 1)
            src_i = mod1(i - sy, nrow)
            field[i, j] *= cis(scale * T(dmap[src_j, src_i]))
        end
    end
    return wf
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

    dm_grid = _dm_influence_grid(dm_cmd, inf, inf_mag, nx_grid, ny_grid, xoff_grid, yoff_grid)

    xdim = min(round(Int, sqrt(2) * nx_grid * dx_inf / dx_surf), n)
    ydim = min(round(Int, sqrt(2) * ny_grid * dx_inf / dx_surf), n)

    a = deg2rad(float(xtilt))
    b = deg2rad(float(ytilt))
    g = deg2rad(float(ztilt))
    grid = Matrix{eltype(dm_grid)}(undef, ydim, xdim)

    if iszero(a) && iszero(b) && iszero(g)
        xstart = (-(xdim ÷ 2) * dx_surf + float(dm_xc) * dx_dm) / dx_inf + xoff_grid0
        ystart = (-(ydim ÷ 2) * dx_surf + float(dm_yc) * dx_dm) / dx_inf + yoff_grid0
        xcoords = range(xstart, step=dx_surf / dx_inf, length=xdim)
        ycoords = range(ystart, step=dx_surf / dx_inf, length=ydim)
        _dm_cubic_conv_transposed_grid!(grid, dm_grid, xcoords, ycoords)
    else
        xax = range(-(xdim ÷ 2) * dx_surf, step=dx_surf, length=xdim)
        yax = range(-(ydim ÷ 2) * dx_surf, step=dx_surf, length=ydim)
        x = repeat(reshape(xax, 1, :), ydim, 1)
        y = repeat(reshape(yax, :, 1), 1, xdim)

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

        _dm_cubic_conv_transposed_coordinate_grid!(grid, dm_grid, xdm, ydm)
    end

    gy, gx = size(grid)
    xmin = n ÷ 2 - xdim ÷ 2 + 1
    ymin = n ÷ 2 - ydim ÷ 2 + 1
    xmax = xmin + gx - 1
    ymax = ymin + gy - 1

    1 <= xmin <= xmax <= n || throw(ArgumentError("PROP_DM: X placement out of bounds"))
    1 <= ymin <= ymax <= n || throw(ArgumentError("PROP_DM: Y placement out of bounds"))
    dmap = if xdim == n && ydim == n && xmin == 1 && ymin == 1
        grid
    else
        out = zeros(eltype(grid), n, n)
        @views out[ymin:ymax, xmin:xmax] .= grid
        out
    end

    if !switch_set(:NO_APPLY; kwargs...)
        # Match accepted D-0036 semantics: the projected DM map must be
        # transposed before centered-map application onto the wavefront grid.
        _dm_apply_projected_phase!(wf, dmap)
    end

    return dmap
end
