using Proper

function telescope_dm(wfo::WaveFront, f_lens::Real, use_errors::Bool, use_dm::Bool)
    obj_map = nothing
    if use_errors
        rms_error = 10e-9
        c_freq = 15.0
        high_power = 3.0
        obj_map = prop_psd_errormap(wfo, rms_error, c_freq, high_power; RMS=true, MAP="obj_map", FILE="telescope_obj.fits")
    end

    prop_lens(wfo, f_lens, "objective")
    prop_propagate(wfo, 2f_lens, "telescope pupil imaging lens")
    prop_lens(wfo, f_lens, "telescope pupil imaging lens")
    prop_propagate(wfo, f_lens, "DM")

    if use_dm && obj_map !== nothing
        nact = 49
        nact_across_pupil = 47
        dm_xc = nact ÷ 2
        dm_yc = nact ÷ 2
        d_beam = 2 * prop_get_beamradius(wfo)
        act_spacing = d_beam / nact_across_pupil
        map_spacing = prop_get_sampling(wfo)

        rot = rotr90(rotr90(obj_map))
        rot = circshift(rot, (1, 1))
        dm_map = prop_magnify(rot, map_spacing / act_spacing, nact; QUICK=true)
        prop_dm(wfo, -dm_map ./ 2, dm_xc, dm_yc, act_spacing; FIT=true)
    end

    prop_propagate(wfo, f_lens, "coronagraph lens")
    return wfo
end
