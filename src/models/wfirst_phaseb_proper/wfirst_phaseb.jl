function _wfirst_phaseb_impl(lambda_m, output_dim0, passvalue; assets=nothing)
    cor_type = String(passget(passvalue, :cor_type, "hlc"))
    use_errors = Int(passget(passvalue, :use_errors, 1))
    use_errors == 0 || throw(ArgumentError("wfirst_phaseb CPU comparison path currently supports only use_errors=0"))
    use_hlc_dm_patterns = Int(passget(passvalue, :use_hlc_dm_patterns, 0))
    use_dm1 = Int(passget(passvalue, :use_dm1, 0))
    use_dm2 = Int(passget(passvalue, :use_dm2, 0))
    dm1_m = passget(passvalue, :dm1_m, nothing)
    dm2_m = passget(passvalue, :dm2_m, nothing)
    dm_sampling_m = Float64(passget(passvalue, :dm_sampling_m, 0.9906e-3))
    dm1_xc_act = Float64(passget(passvalue, :dm1_xc_act, 23.5))
    dm1_yc_act = Float64(passget(passvalue, :dm1_yc_act, 23.5))
    dm1_xtilt_deg = Float64(passget(passvalue, :dm1_xtilt_deg, 0.0))
    dm1_ytilt_deg = Float64(passget(passvalue, :dm1_ytilt_deg, 5.7))
    dm1_ztilt_deg = Float64(passget(passvalue, :dm1_ztilt_deg, 0.0))
    dm2_xc_act = Float64(passget(passvalue, :dm2_xc_act, 23.5))
    dm2_yc_act = Float64(passget(passvalue, :dm2_yc_act, 23.5))
    dm2_xtilt_deg = Float64(passget(passvalue, :dm2_xtilt_deg, 0.0))
    dm2_ytilt_deg = Float64(passget(passvalue, :dm2_ytilt_deg, 5.7))
    dm2_ztilt_deg = Float64(passget(passvalue, :dm2_ztilt_deg, 0.0))
    use_fpm = Int(passget(passvalue, :use_fpm, 1))
    use_pupil_mask = Int(passget(passvalue, :use_pupil_mask, 1))
    use_lyot_stop = Int(passget(passvalue, :use_lyot_stop, 1))
    use_field_stop = Int(passget(passvalue, :use_field_stop, 1))
    final_sampling_m = Float64(passget(passvalue, :final_sampling_m, 0.0))
    final_sampling_lam0 = Float64(passget(passvalue, :final_sampling_lam0, 0.0))
    output_dim = Int(passget(passvalue, :output_dim, output_dim0))
    λm = Float64(lambda_m)
    data_root = String(passget(passvalue, :data_dir, phaseb_default_data_root()))

    diam = 2.3633372
    fl_pri = 2.83459423440 * 1.0013
    d_pri_sec = 2.285150515460035
    fl_sec = -0.653933011 * 1.0004095
    diam_sec = 0.58166
    d_sec_fold1 = 2.993753476654728
    diam_fold1 = 0.09
    d_fold1_m3 = 1.680935841598811
    fl_m3 = 0.430216463069001
    diam_m3 = 0.2
    d_m3_m4 = 0.943514749358944
    fl_m4 = 0.116239114833590
    diam_m4 = 0.07
    d_m4_m5 = 0.429145636743193
    fl_m5 = 0.198821518772608
    diam_m5 = 0.07
    d_m5_fold2 = 0.351125431220770
    diam_fold2 = 0.06
    d_fold2_fsm = 0.365403811661862
    d_fsm_oap1 = 0.354826767220001
    fl_oap1 = 0.503331895563883
    diam_oap1 = 0.06
    d_oap1_focm = 0.768005607094041
    d_focm_oap2 = 0.314483210543378
    fl_oap2 = 0.579156922073536
    diam_oap2 = 0.06
    d_oap2_dm1 = 0.775775726154228
    d_dm1_dm2 = 1.0
    d_dm2_oap3 = 0.394833855161549
    fl_oap3 = 1.217276467668519
    diam_oap3 = 0.06
    d_oap3_fold3 = 0.505329955078121
    diam_fold3 = 0.06
    d_fold3_oap4 = 1.158897671642761
    fl_oap4 = 0.446951159052363
    diam_oap4 = 0.06
    d_oap4_pupilmask = 0.423013568764728
    d_pupilmask_oap5 = 0.408810648253099
    fl_oap5 = 0.548189351937178
    diam_oap5 = 0.06
    d_oap5_fpm = 0.548189083164429
    d_fpm_oap6 = 0.548189083164429
    fl_oap6 = 0.548189083164429
    diam_oap6 = 0.06
    d_oap6_lyotstop = 0.687567667550736
    d_lyotstop_oap7 = 0.401748843470518
    fl_oap7 = 0.708251083480054
    diam_oap7 = 0.06
    d_oap7_fieldstop = 0.708251083480054
    d_fieldstop_oap8 = 0.210985967281651
    fl_oap8 = 0.210985967281651
    diam_oap8 = 0.06
    d_oap8_filter = 0.368452268225530
    diam_filter = 0.01
    d_filter_lens = 0.170799548215162
    fl_lens = 0.246017378417573 + 0.050001306014153
    diam_lens = 0.01
    d_lens_fold4 = 0.246017378417573
    diam_fold4 = 0.02
    d_fold4_image = 0.050001578514650

    cfg = _phaseb_config(cor_type, λm, data_root; compact=false, use_fpm=use_fpm)
    λ0 = cfg.lambda0_m
    pupil_diam_pix = cfg.pupil_diam_pix
    data = cfg.branch == :hlc ? _resolve_assets(passvalue, assets, λm, output_dim) : nothing
    ws = data === nothing ? PhaseBModelWorkspace(output_dim) : data.workspace
    if cfg.branch == :none
        use_fpm = 0
        use_lyot_stop = 0
        use_field_stop = 0
    end
    n_default = cfg.n_default
    n_to_fpm = cfg.n_to_fpm
    n_from_lyotstop = cfg.n_from_lyotstop
    field_stop_radius_lam0 = get(cfg, :field_stop_radius_lam0, 9.0)
    source_x_offset, source_y_offset = _source_offset_lambda_over_d(passvalue, λ0, diam)

    n = n_default
    wf = prop_begin(diam, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    pupil = cfg.branch == :hlc ? data.shared.pupil_1024 : trim(Float64.(prop_fits_read(cfg.pupil_file)), n)
    prop_multiply(wf, pupil)
    prop_define_entrance(wf)
    _apply_source_offset!(wf, pupil_diam_pix, λ0, λm, source_x_offset, source_y_offset)
    prop_lens(wf, fl_pri)

    prop_propagate(wf, d_pri_sec, "secondary")
    prop_lens(wf, fl_sec)

    prop_propagate(wf, d_sec_fold1, "FOLD_1")
    prop_propagate(wf, d_fold1_m3, "M3")
    prop_lens(wf, fl_m3)
    prop_propagate(wf, d_m3_m4, "M4")
    prop_lens(wf, fl_m4)
    prop_propagate(wf, d_m4_m5, "M5")
    prop_lens(wf, fl_m5)
    prop_propagate(wf, d_m5_fold2, "FOLD_2")
    prop_propagate(wf, d_fold2_fsm, "FSM")
    prop_propagate(wf, d_fsm_oap1, "OAP1")
    prop_lens(wf, fl_oap1)
    prop_propagate(wf, d_oap1_focm, "FOCM")
    prop_propagate(wf, d_focm_oap2, "OAP2")
    prop_lens(wf, fl_oap2)

    prop_propagate(wf, d_oap2_dm1, "DM1")
    if use_dm1 != 0
        dm1_m === nothing && throw(ArgumentError("wfirst_phaseb requires dm1_m when use_dm1 != 0"))
        prop_dm(wf, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m; XTILT=dm1_xtilt_deg, YTILT=dm1_ytilt_deg, ZTILT=dm1_ztilt_deg)
    end
    if use_hlc_dm_patterns != 0
        cfg.branch == :hlc || throw(ArgumentError("use_hlc_dm_patterns is only valid for HLC configurations"))
        prop_add_phase(wf, data.shared.dm1wfe_1024)
    end
    prop_propagate(wf, d_dm1_dm2, "DM2")
    if use_dm2 != 0
        dm2_m === nothing && throw(ArgumentError("wfirst_phaseb requires dm2_m when use_dm2 != 0"))
        prop_dm(wf, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m; XTILT=dm2_xtilt_deg, YTILT=dm2_ytilt_deg, ZTILT=dm2_ztilt_deg)
    end
    if use_hlc_dm_patterns != 0
        cfg.branch == :hlc || throw(ArgumentError("use_hlc_dm_patterns is only valid for HLC configurations"))
        prop_add_phase(wf, data.shared.dm2wfe_1024)
    end
    if cfg.branch == :hlc
        prop_multiply(wf, data.shared.dm2mask_1024)
    end

    prop_propagate(wf, d_dm2_oap3, "OAP3")
    prop_lens(wf, fl_oap3)
    prop_propagate(wf, d_oap3_fold3, "FOLD_3")
    prop_propagate(wf, d_fold3_oap4, "OAP4")
    prop_lens(wf, fl_oap4)
    prop_propagate(wf, d_oap4_pupilmask, "PUPIL_MASK")
    if cfg.branch == :spc && use_pupil_mask != 0
        pupil_mask = trim(Float64.(prop_fits_read(cfg.pupil_mask_file)), n)
        prop_multiply(wf, pupil_mask)
    end

    diam_pupil = 2 * prop_get_beamradius(wf)
    field_default = phaseb_field(ws, n_default)
    prop_end!(field_default, wf; noabs=true)
    n = n_to_fpm
    field_to_fpm = phaseb_field(ws, n)
    phaseb_center_copy!(field_to_fpm, field_default)
    wf = prop_begin(diam_pupil, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    phaseb_half_shift!(wf.field, field_to_fpm)

    prop_propagate(wf, d_pupilmask_oap5, "OAP5")
    prop_lens(wf, fl_oap5)
    prop_propagate(wf, d_oap5_fpm, "FPM"; TO_PLANE=true)
    if use_fpm == 1
        if cfg.branch == :hlc
            prop_multiply(wf, data.occulter_2048)
        elseif cfg.branch == :spc
            prop_end!(field_to_fpm, wf; noabs=true)
            phaseb_ffts!(field_to_fpm, phaseb_fft_cache(ws, n), +1)
            field_mft = phaseb_field(ws, cfg.n_mft)
            phaseb_center_copy!(field_mft, field_to_fpm)
            fpm = ComplexF64.(prop_fits_read(cfg.fpm_file))
            nfpm = size(fpm, 2)
            fpm_sampling_lam = cfg.fpm_sampling * cfg.fpm_sampling_lambda_m / λm
            field_spc = mft2(field_mft, fpm_sampling_lam, pupil_diam_pix, nfpm, -1)
            field_spc .*= fpm
            field_to_fpm = mft2(field_spc, fpm_sampling_lam, pupil_diam_pix, n, +1)
            phaseb_ffts!(field_to_fpm, phaseb_fft_cache(ws, n), -1)
            phaseb_half_shift!(wf.field, field_to_fpm)
        end
    end

    prop_propagate(wf, d_fpm_oap6, "OAP6")
    prop_lens(wf, fl_oap6)
    prop_propagate(wf, d_oap6_lyotstop, "LYOT_STOP")

    diam_lyot = 2 * prop_get_beamradius(wf)
    prop_end!(field_to_fpm, wf; noabs=true)
    n = n_from_lyotstop
    field_from_lyot = phaseb_field(ws, n)
    phaseb_center_copy!(field_from_lyot, field_to_fpm)
    wf = prop_begin(diam_lyot, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    phaseb_half_shift!(wf.field, field_from_lyot)

    if use_lyot_stop != 0 && cfg.lyot_stop_file !== nothing
        lyot = cfg.branch == :hlc ? data.shared.lyot_1024 : trim(Float64.(prop_fits_read(cfg.lyot_stop_file)), n)
        prop_multiply(wf, lyot)
    end

    prop_propagate(wf, d_lyotstop_oap7, "OAP7")
    prop_lens(wf, fl_oap7)
    prop_propagate(wf, d_oap7_fieldstop, "FIELD_STOP")
    if use_field_stop != 0 && cfg.branch == :hlc
        sampling_lamD = pupil_diam_pix / n
        stop_radius = field_stop_radius_lam0 / sampling_lamD * (λ0 / λm) * prop_get_sampling(wf)
        prop_circular_aperture(wf, stop_radius, 0.0, 0.0)
    end

    prop_propagate(wf, d_fieldstop_oap8, "OAP8")
    prop_lens(wf, fl_oap8)
    prop_propagate(wf, d_oap8_filter, "filter")
    prop_propagate(wf, d_filter_lens, "LENS")
    prop_lens(wf, fl_lens)
    prop_propagate(wf, d_lens_fold4, "FOLD_4")
    prop_propagate(wf, d_fold4_image, "IMAGE")

    sampling_m = prop_get_sampling(wf)
    prop_end!(field_from_lyot, wf; noabs=true)
    if final_sampling_lam0 != 0 || final_sampling_m != 0
        mag = if final_sampling_m != 0
            sampling_m / final_sampling_m
        else
            (pupil_diam_pix / n) / final_sampling_lam0 * (λm / λ0)
        end
        sampling_m = final_sampling_m != 0 ? final_sampling_m : sampling_m / mag
        ctx = Proper.active_run_context()
        if ctx === nothing
            prop_magnify!(ws.output, field_from_lyot, mag; AMP_CONSERVE=true)
        else
            prop_magnify!(ws.output, field_from_lyot, mag, ctx; AMP_CONSERVE=true)
        end
        return ws.output, sampling_m
    else
        phaseb_center_copy!(ws.output, field_from_lyot)
        return ws.output, sampling_m
    end
end

function wfirst_phaseb(lambda_m, output_dim0, passvalue; assets=nothing)
    return _wfirst_phaseb_impl(lambda_m, output_dim0, passvalue; assets=assets)
end

function wfirst_phaseb(lambda_m, output_dim0; PASSVALUE=nothing, assets=nothing)
    return _wfirst_phaseb_impl(lambda_m, output_dim0, PASSVALUE; assets=assets)
end
