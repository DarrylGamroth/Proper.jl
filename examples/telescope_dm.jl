using Proper

function telescope_dm(wfo::WaveFront, f_lens::Real, use_errors::Bool, use_dm::Bool)
    obj_map = nothing
    if use_errors
        rms_error = 10e-9
        c_freq = 15.0
        high_power = 3.0
        obj_map = prop_psd_errormap(wfo, rms_error, c_freq, high_power; RMS=true, MAP="obj_map", FILE="telescope_obj.fits", no_apply=true)
        prop_add_phase(wfo, (2pi / prop_get_wavelength(wfo)) .* obj_map)
    end

    prop_lens(wfo, f_lens, "objective")
    prop_propagate(wfo, 2f_lens, "telescope pupil imaging lens")
    prop_lens(wfo, f_lens, "telescope pupil imaging lens")
    prop_propagate(wfo, f_lens, "DM")

    if use_dm && obj_map !== nothing
        nact = 49
        nact_across_pupil = 47
        d_beam = 2 * prop_get_beamradius(wfo)
        act_spacing = d_beam / nact_across_pupil
        map_spacing = prop_get_sampling(wfo)

        rot = rotr90(rotr90(obj_map))
        rot = circshift(rot, (1, 1))
        dm_map = prop_magnify(rot, map_spacing / act_spacing, nact)
        prop_dm(wfo, -dm_map ./ 2; mirror=true)
    end

    prop_propagate(wfo, f_lens, "coronagraph lens")
    return wfo
end
