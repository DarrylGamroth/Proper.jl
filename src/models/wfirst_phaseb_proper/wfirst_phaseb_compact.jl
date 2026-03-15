function _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=nothing)
    cor_type = String(passget(passvalue, :cor_type, "hlc"))
    cor_type == "hlc" || throw(ArgumentError("wfirst_phaseb_compact supports only cor_type=hlc for the benchmarked CPU comparison path"))
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
    λ0 = 0.575e-6
    pupil_diam_pix = 309.0
    n_small = 1024
    n_big = 2048
    diam_at_dm1 = 0.0463
    d_dm1_dm2 = 1.0
    source_x_offset, source_y_offset = _source_offset_lambda_over_d(passvalue, λ0, diam_at_dm1)
    data = _resolve_assets(passvalue, assets, λm)
    ws = data.workspace

    n = n_small
    wf = prop_begin(diam_at_dm1, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    prop_multiply(wf, data.shared.pupil_1024)
    prop_define_entrance(wf)
    _apply_source_offset!(wf, pupil_diam_pix, λ0, λm, source_x_offset, source_y_offset)

    if use_dm1 != 0
        dm1_m === nothing && throw(ArgumentError("wfirst_phaseb_compact requires dm1_m when use_dm1 != 0"))
        prop_dm(wf, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m; XTILT=dm1_xtilt_deg, YTILT=dm1_ytilt_deg, ZTILT=dm1_ztilt_deg)
    end

    if use_hlc_dm_patterns != 0
        prop_add_phase(wf, data.shared.dm1wfe_1024)
    end
    prop_propagate(wf, d_dm1_dm2, "DM2")
    if use_dm2 != 0
        dm2_m === nothing && throw(ArgumentError("wfirst_phaseb_compact requires dm2_m when use_dm2 != 0"))
        prop_dm(wf, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m; XTILT=dm2_xtilt_deg, YTILT=dm2_ytilt_deg, ZTILT=dm2_ztilt_deg)
    end
    if use_hlc_dm_patterns != 0
        prop_add_phase(wf, data.shared.dm2wfe_1024)
    end
    prop_multiply(wf, data.shared.dm2mask_1024)
    prop_propagate(wf, -d_dm1_dm2, "back to DM1")

    prop_end!(ws.field_1024, wf; noabs=true)
    phaseb_center_copy!(ws.field_2048, ws.field_1024)
    phaseb_ffts!(ws.field_2048, ws.fft_2048, -1)
    ws.field_2048 .*= data.occulter_2048
    phaseb_ffts!(ws.field_2048, ws.fft_2048, +1)

    phaseb_center_copy!(ws.field_1024, ws.field_2048)
    ws.field_1024 .*= data.shared.lyot_1024
    ws.field_1024 .*= n_small
    phaseb_ffts!(ws.field_1024, ws.fft_1024, -1)
    phaseb_reverse_shift1!(ws.fft_1024.shifted, ws.field_1024)
    copyto!(ws.field_1024, ws.fft_1024.shifted)

    if final_sampling_lam0 != 0
        mag = (pupil_diam_pix / n_small) / final_sampling_lam0 * (λm / λ0)
        ctx = Proper.active_run_context()
        if ctx === nothing
            prop_magnify!(ws.output, ws.field_1024, mag; AMP_CONSERVE=true)
        else
            prop_magnify!(ws.output, ws.field_1024, mag, ctx; AMP_CONSERVE=true)
        end
        return ws.output, 0.0
    else
        phaseb_center_copy!(ws.output, ws.field_1024)
        return ws.output, 0.0
    end
end

function wfirst_phaseb_compact(lambda_m, output_dim0, passvalue; assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=assets)
end

function wfirst_phaseb_compact(lambda_m, output_dim0; PASSVALUE=nothing, assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, PASSVALUE; assets=assets)
end
