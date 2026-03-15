function _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=nothing)
    cor_type = String(passget(passvalue, :cor_type, "hlc"))
    cor_type == "hlc" || throw(ArgumentError("wfirst_phaseb_compact supports only cor_type=hlc for the benchmarked CPU comparison path"))
    use_hlc_dm_patterns = Int(passget(passvalue, :use_hlc_dm_patterns, 0))
    final_sampling_lam0 = Float64(passget(passvalue, :final_sampling_lam0, 0.0))
    output_dim = Int(passget(passvalue, :output_dim, output_dim0))
    λm = Float64(lambda_m)
    λ0 = 0.575e-6
    pupil_diam_pix = 309.0
    n_small = 1024
    n_big = 2048
    diam_at_dm1 = 0.0463
    d_dm1_dm2 = 1.0
    data = _resolve_assets(passvalue, assets, λm)
    ws = data.workspace

    n = n_small
    wf = prop_begin(diam_at_dm1, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    prop_multiply(wf, data.shared.pupil_1024)
    prop_define_entrance(wf)

    if use_hlc_dm_patterns != 0
        prop_add_phase(wf, data.shared.dm1wfe_1024)
    end
    prop_propagate(wf, d_dm1_dm2, "DM2")
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
