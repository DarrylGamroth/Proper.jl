@testset "DM projection active path" begin
    function _dm_projection_dm_command(side::Integer; sparse::Bool)
        n = Int(side)
        cmd = zeros(Float64, n, n)
        center = (n + 1) / 2
        denom = max((n - 1) / 2, 1)
        @inbounds for ix in 1:n
            x = (ix - center) / denom
            for iy in 1:n
                sparse && ((ix + 2iy) % 3 == 0) && continue
                y = (iy - center) / denom
                cmd[iy, ix] = 8e-10 * (1 + 0.07x - 0.04y + 0.02sinpi((iy + ix) / n))
            end
        end
        return cmd
    end

    function _dm_projection_seed_field!(wf)
        @inbounds for j in axes(wf.field, 2)
            for i in axes(wf.field, 1)
                phase = 0.002 * i - 0.0015 * j
                wf.field[i, j] = cis(phase)
            end
        end
        return wf
    end

    function _dm_projection_dm_reference(wf, dm_cmd, dm_xc, dm_yc; nact)
        influence_path = joinpath(@__DIR__, "..", "data", "influence_dm5v2_1.fits")
        inf, hdr = prop_fits_read(influence_path; header=true)
        dx_inf_native = float(hdr["P2PD_M"])
        dx_dm_inf = float(hdr["C2CD_M"])
        inf_mag = round(Int, dx_dm_inf / dx_inf_native)

        n = size(wf.field, 1)
        dx_surf = wf.sampling_m
        dx_dm = 2 * prop_get_beamradius(wf) / float(nact)
        dx_inf = dx_inf_native * dx_dm / dx_dm_inf

        ny_dm, nx_dm = size(dm_cmd)
        margin = 9 * inf_mag
        nx_grid = nx_dm * inf_mag + 2 * margin
        ny_grid = ny_dm * inf_mag + 2 * margin
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
        dm_grid = Proper._fftconvolve_same(dm_grid, inf)

        xdim = min(round(Int, sqrt(2) * nx_grid * dx_inf / dx_surf), n)
        ydim = min(round(Int, sqrt(2) * ny_grid * dx_inf / dx_surf), n)
        xax = range(-(xdim ÷ 2) * dx_surf, step=dx_surf, length=xdim)
        yax = range(-(ydim ÷ 2) * dx_surf, step=dx_surf, length=ydim)
        x = repeat(reshape(xax, 1, :), ydim, 1)
        y = repeat(reshape(yax, :, 1), 1, xdim)
        xdm = (x .+ float(dm_xc) * dx_dm) ./ dx_inf .+ xoff_grid0
        ydm = (y .+ float(dm_yc) * dx_dm) ./ dx_inf .+ yoff_grid0

        # Python passes an F-contiguous `dm_grid.T` view to a row-major C
        # routine. The C routine consumes the underlying `dm_grid` memory, so
        # the equivalent logical Julia oracle samples `dm_grid` directly.
        # Keeping this oracle on the public interpolation path also makes it
        # independent of the private allocation-avoiding DM projection loops.
        grid = prop_cubic_conv(dm_grid, xdm, ydm; grid=false)
        dmap = zeros(eltype(grid), n, n)
        gy, gx = size(grid)
        xmin = n ÷ 2 - xdim ÷ 2 + 1
        ymin = n ÷ 2 - ydim ÷ 2 + 1
        @views dmap[ymin:(ymin + gy - 1), xmin:(xmin + gx - 1)] .= grid
        return dmap
    end

    for (side, sparse, dm_xc, dm_yc) in ((5, true, 0.25, -0.5), (8, false, -0.4, 0.3))
        dm = _dm_projection_dm_command(side; sparse=sparse)

        wf_map = prop_begin(1.0, 0.55e-6, 64)
        fast = prop_dm(wf_map, dm, dm_xc, dm_yc, 0.0; N_ACT_ACROSS_PUPIL=side, NO_APPLY=true)
        ref = _dm_projection_dm_reference(wf_map, dm, dm_xc, dm_yc; nact=side)
        @test fast ≈ ref rtol=5e-13 atol=1e-22

        wf_fast = _dm_projection_seed_field!(prop_begin(1.0, 0.55e-6, 64))
        wf_ref = _dm_projection_seed_field!(prop_begin(1.0, 0.55e-6, 64))
        prop_dm(wf_fast, dm, dm_xc, dm_yc, 0.0; N_ACT_ACROSS_PUPIL=side)
        prop_add_phase(wf_ref, 2 .* ref)
        @test wf_fast.field ≈ wf_ref.field rtol=5e-13 atol=5e-15
    end
end
