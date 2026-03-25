using Proper

@inline function _migration_passget(passvalue, key::Symbol, default=nothing)
    if passvalue === nothing
        return default
    elseif passvalue isa NamedTuple
        return get(passvalue, key, default)
    elseif passvalue isa AbstractDict
        return get(passvalue, key, get(passvalue, String(key), default))
    end
    return default
end

function migration_dm_fits_prescription(λm, n; PASSVALUE=nothing)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.45)
    prop_define_entrance(wf)

    errormap_path = _migration_passget(PASSVALUE, :errormap_path, nothing)
    dm_map = _migration_passget(PASSVALUE, :dm_map, nothing)

    if errormap_path !== nothing
        prop_errormap(wf, errormap_path; WAVEFRONT=true, SAMPLING=prop_get_sampling(wf))
    end

    if dm_map !== nothing
        prop_dm(wf, dm_map)
    end

    prop_lens(wf, 10.0)
    prop_propagate(wf, 10.0)
    return prop_end(wf)
end

function migration_dm_fits_demo(; wavelength_microns::Real=0.55, gridsize::Integer=64)
    mktempdir() do dir
        errormap_path = joinpath(dir, "phase_map.fits")

        errormap = zeros(Float64, gridsize, gridsize)
        errormap[gridsize ÷ 2, gridsize ÷ 2] = 5e-9
        prop_fits_write(errormap_path, errormap)

        dm_map = zeros(Float64, gridsize, gridsize)
        cy = gridsize ÷ 2
        cx = gridsize ÷ 2
        for y in 1:gridsize, x in 1:gridsize
            r2 = (x - cx)^2 + (y - cy)^2
            dm_map[y, x] = 2e-9 * exp(-r2 / max(gridsize, 1))
        end

        passvalue = Dict(
            :errormap_path => errormap_path,
            :dm_map => dm_map,
        )

        return prop_run(migration_dm_fits_prescription, wavelength_microns, gridsize; PASSVALUE=passvalue)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    psf, sampling = migration_dm_fits_demo()
    println("Migration DM/FITS example")
    println("  output size = ", size(psf))
    println("  sampling = ", sampling)
    println("  peak = ", maximum(psf))
end
