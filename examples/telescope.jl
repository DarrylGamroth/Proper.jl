using Proper

function telescope(wfo::WaveFront, f_lens::Real, use_errors::Bool, use_dm::Bool=false)
    if use_errors
        rms_error = 10e-9
        c_freq = 15.0
        high_power = 3.0
        prop_psd_errormap(wfo, rms_error, c_freq, high_power; RMS=true, MAP="obj_map", FILE="telescope_obj.fits")
    end

    prop_lens(wfo, f_lens, "objective")
    prop_propagate(wfo, 2f_lens, "telescope pupil imaging lens")
    prop_lens(wfo, f_lens, "telescope pupil imaging lens")
    prop_propagate(wfo, f_lens, "DM")
    prop_propagate(wfo, f_lens, "coronagraph lens")
    return wfo
end
