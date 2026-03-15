function _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=nothing)
    cor_type = String(passget(passvalue, :cor_type, "hlc"))
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
    final_sampling_lam0 = Float64(passget(passvalue, :final_sampling_lam0, 0.0))
    output_dim = Int(passget(passvalue, :output_dim, output_dim0))
    λm = Float64(lambda_m)
    diam_at_dm1 = 0.0463
    d_dm1_dm2 = 1.0
    data_root = String(passget(passvalue, :data_dir, phaseb_default_data_root()))
    cfg = _phaseb_config(cor_type, λm, data_root; compact=true)
    λ0 = cfg.lambda0_m
    pupil_diam_pix = cfg.pupil_diam_pix
    n_small = cfg.n_small
    n_big = cfg.n_big
    source_x_offset, source_y_offset = _source_offset_lambda_over_d(passvalue, λ0, diam_at_dm1)
    data = cfg.branch == :hlc ? _resolve_assets(passvalue, assets, λm) : nothing
    ws = data === nothing ? PhaseBModelWorkspace(output_dim) : data.workspace
    field_small = phaseb_field(ws, n_small)
    fft_small = phaseb_fft_cache(ws, n_small)

    n = n_small
    wf = prop_begin(diam_at_dm1, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    pupil = cfg.branch == :hlc ? data.shared.pupil_1024 : trim(Float64.(prop_fits_read(cfg.pupil_file)), n)
    prop_multiply(wf, pupil)
    prop_define_entrance(wf)
    _apply_source_offset!(wf, pupil_diam_pix, λ0, λm, source_x_offset, source_y_offset)

    if use_dm1 != 0
        dm1_m === nothing && throw(ArgumentError("wfirst_phaseb_compact requires dm1_m when use_dm1 != 0"))
        prop_dm(wf, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m; XTILT=dm1_xtilt_deg, YTILT=dm1_ytilt_deg, ZTILT=dm1_ztilt_deg)
    end

    if use_hlc_dm_patterns != 0
        cfg.branch == :hlc || throw(ArgumentError("use_hlc_dm_patterns is only valid for HLC compact configurations"))
        prop_add_phase(wf, data.shared.dm1wfe_1024)
    end
    prop_propagate(wf, d_dm1_dm2, "DM2")
    if use_dm2 != 0
        dm2_m === nothing && throw(ArgumentError("wfirst_phaseb_compact requires dm2_m when use_dm2 != 0"))
        prop_dm(wf, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m; XTILT=dm2_xtilt_deg, YTILT=dm2_ytilt_deg, ZTILT=dm2_ztilt_deg)
    end
    if use_hlc_dm_patterns != 0
        cfg.branch == :hlc || throw(ArgumentError("use_hlc_dm_patterns is only valid for HLC compact configurations"))
        prop_add_phase(wf, data.shared.dm2wfe_1024)
    end
    if cfg.branch == :hlc
        prop_multiply(wf, data.shared.dm2mask_1024)
    end
    prop_propagate(wf, -d_dm1_dm2, "back to DM1")

    prop_end!(field_small, wf; noabs=true)
    if cfg.branch == :spc
        pupil_mask = trim(Float64.(prop_fits_read(cfg.pupil_mask_file)), n_small)
        field_small .*= pupil_mask
    end
    if cfg.branch == :hlc
        field_big = phaseb_field(ws, n_big)
        fft_big = phaseb_fft_cache(ws, n_big)
        phaseb_center_copy!(field_big, field_small)
        phaseb_ffts!(field_big, fft_big, -1)
        field_big .*= data.occulter_2048
        phaseb_ffts!(field_big, fft_big, +1)
        phaseb_center_copy!(field_small, field_big)
    elseif cfg.branch == :spc
        fpm = ComplexF64.(prop_fits_read(cfg.fpm_file))
        nfpm = size(fpm, 2)
        fpm_sampling_lam = cfg.fpm_sampling * cfg.fpm_sampling_lambda_m / λm
        field_big = phaseb_field(ws, n_big)
        phaseb_center_copy!(field_big, field_small)
        spc_field = mft2(field_big, fpm_sampling_lam, pupil_diam_pix, nfpm, -1)
        spc_field .*= fpm
        spc_pupil_diam_pix = pupil_diam_pix / 2
        field_big = mft2(spc_field, fpm_sampling_lam, spc_pupil_diam_pix, Int(round(spc_pupil_diam_pix)), +1)
        pupil_diam_pix = spc_pupil_diam_pix
        phaseb_center_copy!(field_small, field_big)
    end

    if cfg.lyot_stop_file !== nothing
        lyot = cfg.branch == :hlc ? data.shared.lyot_1024 : trim(Float64.(prop_fits_read(cfg.lyot_stop_file)), n_small)
        field_small .*= lyot
    end
    field_small .*= n_small
    phaseb_ffts!(field_small, fft_small, -1)
    phaseb_reverse_shift1!(fft_small.shifted, field_small)
    copyto!(field_small, fft_small.shifted)

    if final_sampling_lam0 != 0
        mag = (pupil_diam_pix / n_small) / final_sampling_lam0 * (λm / λ0)
        ctx = Proper.active_run_context()
        if ctx === nothing
            prop_magnify!(ws.output, field_small, mag; AMP_CONSERVE=true)
        else
            prop_magnify!(ws.output, field_small, mag, ctx; AMP_CONSERVE=true)
        end
        return ws.output, 0.0
    else
        phaseb_center_copy!(ws.output, field_small)
        return ws.output, 0.0
    end
end

function wfirst_phaseb_compact(lambda_m, output_dim0, passvalue; assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=assets)
end

function wfirst_phaseb_compact(lambda_m, output_dim0; PASSVALUE=nothing, assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, PASSVALUE; assets=assets)
end
