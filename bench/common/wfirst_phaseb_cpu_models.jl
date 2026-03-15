module WFIRSTPhaseBCPUModels

using FFTW
using LinearAlgebra
using Proper

export phaseb_case_definitions, phaseb_default_data_root, load_phaseb_hlc_assets
export prepare_phaseb_models, run_phaseb_case
export wfirst_phaseb_compact, wfirst_phaseb

struct PhaseBHLCAssets
    data_root::String
    pupil::Matrix{Float64}
    lyot::Matrix{Float64}
    dm1wfe::Matrix{Float64}
    dm2wfe::Matrix{Float64}
    dm2mask::Matrix{Float64}
    occulters::Dict{Float64, Matrix{ComplexF64}}
end

struct PhaseBHLCPreparedSharedAssets
    data_root::String
    pupil_1024::Matrix{Float64}
    lyot_1024::Matrix{Float64}
    dm1wfe_1024::Matrix{Float64}
    dm2wfe_1024::Matrix{Float64}
    dm2mask_1024::Matrix{Float64}
    occulters_2048::Dict{Float64, Matrix{ComplexF64}}
end

const PhaseBFFTCache = Proper.CenteredFFTCache{Float64}

mutable struct PhaseBModelWorkspace
    field_1024::Matrix{ComplexF64}
    field_2048::Matrix{ComplexF64}
    output::Matrix{ComplexF64}
    fft_1024::PhaseBFFTCache
    fft_2048::PhaseBFFTCache
end

struct PhaseBPreparedAssets
    shared::PhaseBHLCPreparedSharedAssets
    occulter_2048::Matrix{ComplexF64}
    workspace::PhaseBModelWorkspace
end

const OLD_LAM_OCCS = [
    "5.4625e-07", "5.49444444444e-07", "5.52638888889e-07", "5.534375e-07", "5.55833333333e-07", "5.59027777778e-07",
    "5.60625e-07", "5.62222222222e-07", "5.65416666667e-07", "5.678125e-07", "5.68611111111e-07", "5.71805555556e-07",
    "5.75e-07", "5.78194444444e-07", "5.81388888889e-07", "5.821875e-07", "5.84583333333e-07", "5.87777777778e-07",
    "5.89375e-07", "5.90972222222e-07", "5.94166666667e-07", "5.965625e-07", "5.97361111111e-07", "6.00555555556e-07", "6.0375e-07",
]
const OLD_LAM_OCCS_M = parse.(Float64, OLD_LAM_OCCS)

@inline passget(::Nothing, key::Symbol, default) = default
@inline function passget(pass::AbstractDict, key::Symbol, default)
    haskey(pass, key) && return pass[key]
    skey = String(key)
    haskey(pass, skey) && return pass[skey]
    return default
end
@inline function passget(pass::NamedTuple, key::Symbol, default)
    return hasproperty(pass, key) ? getproperty(pass, key) : default
end

@inline function phaseb_default_data_root()
    return abspath(get(ENV, "WFIRST_PHASEB_DATA_ROOT", joinpath(pwd(), ".cache", "wfirst", "data_phaseb_from_roman_preflight")))
end

@inline function phaseb_trim(input_image::AbstractMatrix, output_dim::Integer)
    input_dim = size(input_image, 2)
    if input_dim == output_dim
        return input_image
    elseif output_dim < input_dim
        x1 = input_dim ÷ 2 - output_dim ÷ 2 + 1
        x2 = x1 + output_dim - 1
        return copy(@view input_image[x1:x2, x1:x2])
    end

    output_image = zeros(eltype(input_image), output_dim, output_dim)
    x1 = output_dim ÷ 2 - input_dim ÷ 2 + 1
    x2 = x1 + input_dim - 1
    @views output_image[x1:x2, x1:x2] .= input_image
    return output_image
end

@inline phaseb_center_copy!(out::AbstractMatrix{T}, input_image::AbstractMatrix) where {T} = Proper.copy_centered!(out, input_image)

@inline function phaseb_prepare_static(input_image::AbstractMatrix, output_dim::Integer, ::Type{T}) where {T}
    out = Matrix{T}(undef, output_dim, output_dim)
    phaseb_center_copy!(out, input_image)
    return out
end

@inline function phaseb_prepare_complex(input_image::AbstractMatrix, output_dim::Integer)
    out = Matrix{ComplexF64}(undef, output_dim, output_dim)
    phaseb_center_copy!(out, input_image)
    return out
end

@inline phaseb_half_shift!(out::AbstractMatrix{T}, input::AbstractMatrix) where {T} = Proper.half_shift_copy!(out, input)
@inline phaseb_reverse_shift1!(out::AbstractMatrix{ComplexF64}, input::AbstractMatrix{<:Complex}) = Proper.reverse_shift1!(out, input)

function PhaseBFFTCache(n::Integer; flags=FFTW.MEASURE)
    return Proper.CenteredFFTCache(Float64, n; flags=flags)
end

function PhaseBModelWorkspace(output_dim::Integer)
    return PhaseBModelWorkspace(
        Matrix{ComplexF64}(undef, 1024, 1024),
        Matrix{ComplexF64}(undef, 2048, 2048),
        Matrix{ComplexF64}(undef, output_dim, output_dim),
        PhaseBFFTCache(1024),
        PhaseBFFTCache(2048),
    )
end

function phaseb_ffts!(field::Matrix{ComplexF64}, cache::PhaseBFFTCache, direction::Integer)
    return Proper.centered_fft!(field, cache, direction)
end

function phaseb_ffts(wavefront::AbstractMatrix, direction::Integer)
    arr = ComplexF64.(wavefront)
    n = size(arr, 1)
    shifted = circshift(arr, (-n ÷ 2, -n ÷ 2))
    transformed = if direction == -1
        fft(shifted) ./ length(shifted)
    else
        ifft(shifted) .* length(shifted)
    end
    return circshift(transformed, (n ÷ 2, n ÷ 2))
end

function phaseb_mft2(field_in::AbstractMatrix, dout::Real, D::Real, nout::Integer, direction::Integer; xoffset::Real=0, yoffset::Real=0, xc::Real=0, yc::Real=0)
    nfield_in = size(field_in, 2)
    x = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2) .- xc)
    y = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2) .- yc)
    u = (collect(0:(nout - 1)) .- (nout ÷ 2) .- xoffset / dout) .* (dout / D)
    v = (collect(0:(nout - 1)) .- (nout ÷ 2) .- yoffset / dout) .* (dout / D)
    xu = x * transpose(u)
    yv = y * transpose(v)
    expxu = (dout / D) .* exp.(direction * 2π * im .* xu)
    expyv = transpose(exp.(direction * 2π * im .* yv))
    return expyv * field_in * expxu
end

function _requested_occ_label(λm::Real)
    idx = argmin(abs.(OLD_LAM_OCCS_M .- Float64(λm)))
    return OLD_LAM_OCCS[idx]
end

function load_phaseb_hlc_assets(data_root::AbstractString, wavelengths_m::AbstractVector{<:Real})
    hlc_dir = joinpath(data_root, "hlc_20190210")
    pupil = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_pupil_rotated.fits")))
    lyot = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_lyot.fits")))
    dm1wfe = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_dm1wfe.fits")))
    dm2wfe = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_dm2wfe.fits")))
    dm2mask = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_dm2mask.fits")))
    occulters = Dict{Float64, Matrix{ComplexF64}}()
    for λm0 in wavelengths_m
        λm = Float64(λm0)
        label = _requested_occ_label(λm)
        real_part = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_occ_lam$(label)theta6.69polp_real_rotated.fits")))
        imag_part = Float64.(prop_fits_read(joinpath(hlc_dir, "run461_occ_lam$(label)theta6.69polp_imag_rotated.fits")))
        occulters[λm] = ComplexF64.(real_part .+ im .* imag_part)
    end
    return PhaseBHLCAssets(String(data_root), pupil, lyot, dm1wfe, dm2wfe, dm2mask, occulters)
end

function prepare_phaseb_shared_assets(data_root::AbstractString, wavelengths_m::AbstractVector{<:Real})
    raw = load_phaseb_hlc_assets(data_root, wavelengths_m)
    occulters_2048 = Dict{Float64, Matrix{ComplexF64}}()
    for (λm, occ) in raw.occulters
        occulters_2048[Float64(λm)] = phaseb_prepare_complex(occ, 2048)
    end
    return PhaseBHLCPreparedSharedAssets(
        raw.data_root,
        phaseb_prepare_static(raw.pupil, 1024, Float64),
        phaseb_prepare_static(raw.lyot, 1024, Float64),
        phaseb_prepare_static(raw.dm1wfe, 1024, Float64),
        phaseb_prepare_static(raw.dm2wfe, 1024, Float64),
        phaseb_prepare_static(raw.dm2mask, 1024, Float64),
        occulters_2048,
    )
end

@inline function _resolve_assets(passvalue, assets, λm::Real)
    if assets !== nothing
        return assets
    end
    data_root = String(passget(passvalue, :data_dir, phaseb_default_data_root()))
    shared = prepare_phaseb_shared_assets(data_root, [λm])
    return PhaseBPreparedAssets(shared, _nearest_occulter(shared, λm), PhaseBModelWorkspace(128))
end

@inline function _nearest_occulter(assets::PhaseBHLCAssets, λm::Real)
    keys_vec = collect(keys(assets.occulters))
    key = keys_vec[argmin(abs.(keys_vec .- Float64(λm)))]
    return assets.occulters[key]
end

@inline function _nearest_occulter(assets::PhaseBHLCPreparedSharedAssets, λm::Real)
    keys_vec = collect(keys(assets.occulters_2048))
    key = keys_vec[argmin(abs.(keys_vec .- Float64(λm)))]
    return assets.occulters_2048[key]
end

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

function wfirst_phaseb_compact(lambda_m, output_dim0, passvalue=nothing; assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, passvalue; assets=assets)
end

function wfirst_phaseb_compact(lambda_m, output_dim0; PASSVALUE=nothing, assets=nothing)
    return _wfirst_phaseb_compact_impl(lambda_m, output_dim0, PASSVALUE; assets=assets)
end

function _wfirst_phaseb_impl(lambda_m, output_dim0, passvalue; assets=nothing)
    cor_type = String(passget(passvalue, :cor_type, "hlc"))
    cor_type == "hlc" || throw(ArgumentError("wfirst_phaseb supports only cor_type=hlc for the benchmarked CPU comparison path"))
    use_errors = Int(passget(passvalue, :use_errors, 1))
    use_errors == 0 || throw(ArgumentError("wfirst_phaseb CPU comparison path currently supports only use_errors=0"))
    use_hlc_dm_patterns = Int(passget(passvalue, :use_hlc_dm_patterns, 0))
    final_sampling_lam0 = Float64(passget(passvalue, :final_sampling_lam0, 0.0))
    output_dim = Int(passget(passvalue, :output_dim, output_dim0))
    λm = Float64(lambda_m)
    λ0 = 0.575e-6
    data = _resolve_assets(passvalue, assets, λm)
    ws = data.workspace

    pupil_diam_pix = 309.0
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

    n_default = 1024
    n_to_fpm = 2048
    n_from_lyotstop = 1024
    field_stop_radius_lam0 = 9.0
    use_fpm = 1
    use_lyot_stop = 1
    use_field_stop = 1

    n = n_default
    wf = prop_begin(diam, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    prop_multiply(wf, data.shared.pupil_1024)
    prop_define_entrance(wf)
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
    if use_hlc_dm_patterns != 0
        prop_add_phase(wf, data.shared.dm1wfe_1024)
    end
    prop_propagate(wf, d_dm1_dm2, "DM2")
    if use_hlc_dm_patterns != 0
        prop_add_phase(wf, data.shared.dm2wfe_1024)
    end
    prop_multiply(wf, data.shared.dm2mask_1024)

    prop_propagate(wf, d_dm2_oap3, "OAP3")
    prop_lens(wf, fl_oap3)
    prop_propagate(wf, d_oap3_fold3, "FOLD_3")
    prop_propagate(wf, d_fold3_oap4, "OAP4")
    prop_lens(wf, fl_oap4)
    prop_propagate(wf, d_oap4_pupilmask, "PUPIL_MASK")

    diam_pupil = 2 * prop_get_beamradius(wf)
    prop_end!(ws.field_1024, wf; noabs=true)
    n = n_to_fpm
    phaseb_center_copy!(ws.field_2048, ws.field_1024)
    wf = prop_begin(diam_pupil, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    phaseb_half_shift!(wf.field, ws.field_2048)

    prop_propagate(wf, d_pupilmask_oap5, "OAP5")
    prop_lens(wf, fl_oap5)
    prop_propagate(wf, d_oap5_fpm, "FPM"; TO_PLANE=true)
    if use_fpm == 1
        prop_multiply(wf, data.occulter_2048)
    end

    prop_propagate(wf, d_fpm_oap6, "OAP6")
    prop_lens(wf, fl_oap6)
    prop_propagate(wf, d_oap6_lyotstop, "LYOT_STOP")

    diam_lyot = 2 * prop_get_beamradius(wf)
    prop_end!(ws.field_2048, wf; noabs=true)
    n = n_from_lyotstop
    phaseb_center_copy!(ws.field_1024, ws.field_2048)
    wf = prop_begin(diam_lyot, λm, n; beam_diam_fraction=pupil_diam_pix / n)
    phaseb_half_shift!(wf.field, ws.field_1024)

    if use_lyot_stop != 0
        prop_multiply(wf, data.shared.lyot_1024)
    end

    prop_propagate(wf, d_lyotstop_oap7, "OAP7")
    prop_lens(wf, fl_oap7)
    prop_propagate(wf, d_oap7_fieldstop, "FIELD_STOP")
    if use_field_stop != 0
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
    prop_end!(ws.field_1024, wf; noabs=true)
    if final_sampling_lam0 != 0
        mag = (pupil_diam_pix / n) / final_sampling_lam0 * (λm / λ0)
        sampling_m /= mag
        ctx = Proper.active_run_context()
        if ctx === nothing
            prop_magnify!(ws.output, ws.field_1024, mag; AMP_CONSERVE=true)
        else
            prop_magnify!(ws.output, ws.field_1024, mag, ctx; AMP_CONSERVE=true)
        end
        return ws.output, sampling_m
    else
        phaseb_center_copy!(ws.output, ws.field_1024)
        return ws.output, sampling_m
    end
end

function wfirst_phaseb(lambda_m, output_dim0, passvalue=nothing; assets=nothing)
    return _wfirst_phaseb_impl(lambda_m, output_dim0, passvalue; assets=assets)
end

function wfirst_phaseb(lambda_m, output_dim0; PASSVALUE=nothing, assets=nothing)
    return _wfirst_phaseb_impl(lambda_m, output_dim0, PASSVALUE; assets=assets)
end

function phaseb_case_definitions()
    lam0_um = 0.575
    band = 0.1
    wavelengths_um = collect(range(lam0_um * (1 - band / 2), lam0_um * (1 + band / 2); length=3))
    wavelengths_m = wavelengths_um .* 1.0e-6
    return Dict(
        "compact_hlc" => (
            func=wfirst_phaseb_compact,
            output_dim=128,
            wavelengths_um=wavelengths_um,
            wavelengths_m=wavelengths_m,
            passvalue=Dict(
                "cor_type" => "hlc",
                "use_hlc_dm_patterns" => 1,
                "final_sampling_lam0" => 0.1,
            ),
            description="Compact HLC model with default DM patterns over 10% band",
        ),
        "full_hlc" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=wavelengths_um,
            wavelengths_m=wavelengths_m,
            passvalue=Dict(
                "cor_type" => "hlc",
                "use_hlc_dm_patterns" => 1,
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
            ),
            description="Full HLC model with default DM patterns over 10% band",
        ),
    )
end

function prepare_phaseb_models(case::NamedTuple; data_root::AbstractString=phaseb_default_data_root())
    shared = prepare_phaseb_shared_assets(data_root, case.wavelengths_m)
    passvalue = merge(case.passvalue, Dict("data_dir" => String(data_root)))
    models = [
        prepare_model(
            Symbol("wfirst_" * string(i)),
            case.func,
            λum,
            case.output_dim;
            PASSVALUE=passvalue,
            assets=prepare_asset_pool(() -> PhaseBPreparedAssets(shared, _nearest_occulter(shared, λm), PhaseBModelWorkspace(case.output_dim)); pool_size=1),
            pool_size=1,
        )
        for (i, (λum, λm)) in enumerate(zip(case.wavelengths_um, case.wavelengths_m))
    ]
    return models, shared
end

function run_phaseb_case(models::Vector{<:PreparedModel}; threaded::Bool=true)
    n = length(models)
    outputs = Vector{Matrix{ComplexF64}}(undef, n)
    samplings = Vector{Float64}(undef, n)
    if threaded && Threads.nthreads() > 1 && n > 1
        Threads.@threads for i in 1:n
            out, sampling = prop_run(models[i])
            outputs[i] = ComplexF64.(out)
            samplings[i] = Float64(sampling)
        end
    else
        for i in 1:n
            out, sampling = prop_run(models[i])
            outputs[i] = ComplexF64.(out)
            samplings[i] = Float64(sampling)
        end
    end

    sy, sx = size(outputs[1])
    stack = Array{ComplexF64}(undef, n, sy, sx)
    @inbounds for i in 1:n
        stack[i, :, :] = outputs[i]
    end
    return stack, samplings
end

end
