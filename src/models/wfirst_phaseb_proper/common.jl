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
    fields::Dict{Int, Matrix{ComplexF64}}
    ffts::Dict{Int, PhaseBFFTCache}
    output::Matrix{ComplexF64}
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
const ERKIN_LAM_OCCS = [
    "5.4625e-07", "5.4944e-07", "5.5264e-07", "5.5583e-07", "5.5903e-07", "5.6222e-07", "5.6542e-07",
    "5.6861e-07", "5.7181e-07", "5.75e-07", "5.7819e-07", "5.8139e-07", "5.8458e-07", "5.8778e-07",
    "5.9097e-07", "5.9417e-07", "5.9736e-07", "6.0056e-07", "6.0375e-07",
]
const ERKIN_LAM_OCCS_M = parse.(Float64, ERKIN_LAM_OCCS)

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
        Dict(
            1024 => Matrix{ComplexF64}(undef, 1024, 1024),
            2048 => Matrix{ComplexF64}(undef, 2048, 2048),
        ),
        Dict(
            1024 => PhaseBFFTCache(1024),
            2048 => PhaseBFFTCache(2048),
        ),
        Matrix{ComplexF64}(undef, output_dim, output_dim),
    )
end

@inline function phaseb_field(ws::PhaseBModelWorkspace, n::Integer)
    n > 0 || throw(ArgumentError("grid size must be positive"))
    return get!(ws.fields, Int(n)) do
        Matrix{ComplexF64}(undef, Int(n), Int(n))
    end
end

@inline function phaseb_fft_cache(ws::PhaseBModelWorkspace, n::Integer)
    n > 0 || throw(ArgumentError("grid size must be positive"))
    return get!(ws.ffts, Int(n)) do
        PhaseBFFTCache(Int(n))
    end
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
    expxu = (dout / D) .* exp.(-direction * 2π * im .* xu)
    expyv = transpose(exp.(-direction * 2π * im .* yv))
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

@inline function _requested_occ_label(labels::AbstractVector{<:AbstractString}, labels_m::AbstractVector{<:Real}, λm::Real)
    idx = argmin(abs.(labels_m .- Float64(λm)))
    return labels[idx]
end

@inline function _phaseb_config(cor_type::AbstractString, λm::Real, data_root::AbstractString; compact::Bool=false, use_fpm::Integer=1)
    if cor_type == "hlc"
        hlc_dir = joinpath(data_root, "hlc_20190210")
        label = _requested_occ_label(OLD_LAM_OCCS, OLD_LAM_OCCS_M, λm)
        return (
            branch=:hlc,
            data_dir=hlc_dir,
            pupil_diam_pix=309.0,
            lambda0_m=0.575e-6,
            pupil_file=joinpath(hlc_dir, "run461_pupil_rotated.fits"),
            lyot_stop_file=joinpath(hlc_dir, "run461_lyot.fits"),
            occulter_real_file=joinpath(hlc_dir, "run461_occ_lam$(label)theta6.69polp_real_rotated.fits"),
            occulter_imag_file=joinpath(hlc_dir, "run461_occ_lam$(label)theta6.69polp_imag_rotated.fits"),
            n_small=1024,
            n_big=compact ? 2048 : (use_fpm != 0 ? 2048 : 1024),
            n_default=1024,
            n_to_fpm=use_fpm != 0 ? 2048 : 1024,
            n_from_lyotstop=1024,
            field_stop_radius_lam0=9.0,
        )
    elseif cor_type == "hlc_erkin"
        hlc_dir = joinpath(data_root, "hlc_20190206_v3")
        label = _requested_occ_label(ERKIN_LAM_OCCS, ERKIN_LAM_OCCS_M, λm)
        return (
            branch=:hlc,
            data_dir=hlc_dir,
            pupil_diam_pix=compact ? 310.0 : 310.0,
            lambda0_m=0.575e-6,
            pupil_file=joinpath(hlc_dir, compact ? "dsn17d_run2_pup310_fpm2048_pupil.fits" : "dsn17d_run2_pup310_fpm2048_pupil.fits"),
            lyot_stop_file=joinpath(hlc_dir, compact ? "dsn17d_run2_pup310_fpm2048_lyot.fits" : "dsn17d_run2_pup310_fpm2048_lyot.fits"),
            occulter_real_file=joinpath(hlc_dir, compact ? "dsn17d_run2_pup310_fpm2048_occ_lam$(label)theta6.69pols_real.fits" : "dsn17d_run2_pup310_fpm2048_occ_lam$(label)theta6.69pols_real_rotated.fits"),
            occulter_imag_file=joinpath(hlc_dir, compact ? "dsn17d_run2_pup310_fpm2048_occ_lam$(label)theta6.69pols_imag.fits" : "dsn17d_run2_pup310_fpm2048_occ_lam$(label)theta6.69pols_imag_rotated.fits"),
            n_small=1024,
            n_big=2048,
            n_default=1024,
            n_to_fpm=use_fpm != 0 ? 2048 : 1024,
            n_from_lyotstop=1024,
            field_stop_radius_lam0=9.0,
        )
    elseif cor_type in ("spc-ifs_short", "spc-ifs_long", "spc-spec_short", "spc-spec_long")
        spc_dir = joinpath(data_root, "spc_20190130")
        is_short = cor_type in ("spc-ifs_short", "spc-spec_short")
        return (
            branch=:spc,
            data_dir=spc_dir,
            pupil_diam_pix=1000.0,
            lambda0_m=is_short ? 0.66e-6 : 0.73e-6,
            pupil_file=joinpath(spc_dir, "pupil_SPC-20190130_rotated.fits"),
            pupil_mask_file=joinpath(spc_dir, compact ? "SPM_SPC-20190130_rotated.fits" : "SPM_SPC-20190130.fits"),
            fpm_file=joinpath(spc_dir, "fpm_0.05lamdivD.fits"),
            lyot_stop_file=joinpath(spc_dir, compact ? "lyotstop_0.5mag.fits" : "LS_SPC-20190130.fits"),
            fpm_sampling=0.05,
            fpm_sampling_lambda_m=is_short ? 0.66e-6 : 0.73e-6,
            n_small=2048,
            n_big=compact ? 1400 : 2048,
            n_default=2048,
            n_to_fpm=2048,
            n_mft=1400,
            n_from_lyotstop=compact ? 2048 : 4096,
        )
    elseif cor_type == "spc-wide"
        spc_dir = joinpath(data_root, "spc_20181220")
        return (
            branch=:spc,
            data_dir=spc_dir,
            pupil_diam_pix=1000.0,
            lambda0_m=0.825e-6,
            pupil_file=joinpath(spc_dir, "pupil_SPC-20181220_1k_rotated.fits"),
            pupil_mask_file=joinpath(spc_dir, compact ? "SPM_SPC-20181220_1000_rounded9_gray_rotated.fits" : "SPM_SPC-20181220_1000_rounded9_gray.fits"),
            fpm_file=joinpath(spc_dir, "fpm_0.05lamdivD.fits"),
            lyot_stop_file=joinpath(spc_dir, compact ? "LS_half_symm_CGI180718_Str3.20pct_38D91_N500_pixel.fits" : "LS_SPC-20181220_1k.fits"),
            fpm_sampling=0.05,
            fpm_sampling_lambda_m=0.825e-6,
            n_small=2048,
            n_big=1400,
            n_default=2048,
            n_to_fpm=2048,
            n_mft=1400,
            n_from_lyotstop=compact ? 2048 : 4096,
        )
    elseif cor_type == "none"
        hlc_dir = joinpath(data_root, "hlc_20190210")
        return (
            branch=:none,
            data_dir=hlc_dir,
            pupil_diam_pix=309.0,
            lambda0_m=0.575e-6,
            pupil_file=joinpath(hlc_dir, "run461_pupil_rotated.fits"),
            lyot_stop_file=nothing,
            n_small=1024,
            n_big=1024,
            n_default=1024,
            n_to_fpm=1024,
            n_from_lyotstop=1024,
        )
    end
    throw(ArgumentError("Unsupported cor_type: $(cor_type)"))
end

@inline function _source_offset_lambda_over_d(passvalue, lambda0_m::Real, diam_m::Real)
    source_x_offset = Float64(passget(passvalue, :source_x_offset, 0.0))
    source_y_offset = Float64(passget(passvalue, :source_y_offset, 0.0))
    mas_per_lamD = lambda0_m * 360.0 * 3600.0 / (2π * diam_m) * 1000
    source_x_offset_mas = passget(passvalue, :source_x_offset_mas, nothing)
    source_y_offset_mas = passget(passvalue, :source_y_offset_mas, nothing)
    if source_x_offset_mas !== nothing
        source_x_offset = Float64(source_x_offset_mas) / mas_per_lamD
    end
    if source_y_offset_mas !== nothing
        source_y_offset = Float64(source_y_offset_mas) / mas_per_lamD
    end
    return source_x_offset, source_y_offset
end

function _apply_source_offset!(wf, pupil_diam_pix::Real, lambda0_m::Real, lambda_m::Real, source_x_offset::Real, source_y_offset::Real)
    if source_x_offset == 0 && source_y_offset == 0
        return wf
    end

    xtilt_lam = -Float64(source_x_offset) * lambda0_m / lambda_m
    ytilt_lam = -Float64(source_y_offset) * lambda0_m / lambda_m
    n = size(wf.field, 1)
    coords = (collect(0:(n - 1)) .- (n ÷ 2)) ./ (float(pupil_diam_pix) / 2.0)
    x = repeat(reshape(coords, 1, :), n, 1)
    y = repeat(reshape(coords, :, 1), 1, n)
    wf.field .*= cis.(π .* (xtilt_lam .* x .+ ytilt_lam .* y))
    return wf
end


function phaseb_case_definitions()
    hlc_lam0_um = 0.575
    hlc_band = 0.1
    hlc_wavelengths_um = collect(range(hlc_lam0_um * (1 - hlc_band / 2), hlc_lam0_um * (1 + hlc_band / 2); length=3))
    hlc_wavelengths_m = hlc_wavelengths_um .* 1.0e-6
    spec_lam0_um = 0.73
    spec_band = 0.15
    spec_wavelengths_um = collect(range(spec_lam0_um * (1 - spec_band / 2), spec_lam0_um * (1 + spec_band / 2); length=3))
    spec_wavelengths_m = spec_wavelengths_um .* 1.0e-6
    wide_lam0_um = 0.825
    wide_band = 0.1
    wide_wavelengths_um = collect(range(wide_lam0_um * (1 - wide_band / 2), wide_lam0_um * (1 + wide_band / 2); length=3))
    wide_wavelengths_m = wide_wavelengths_um .* 1.0e-6
    return Dict(
        "compact_hlc" => (
            func=wfirst_phaseb_compact,
            output_dim=128,
            wavelengths_um=hlc_wavelengths_um,
            wavelengths_m=hlc_wavelengths_m,
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
            wavelengths_um=hlc_wavelengths_um,
            wavelengths_m=hlc_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "hlc",
                "use_hlc_dm_patterns" => 1,
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
            ),
            description="Full HLC model with default DM patterns over 10% band",
        ),
        "compact_spc_spec_long" => (
            func=wfirst_phaseb_compact,
            output_dim=128,
            wavelengths_um=spec_wavelengths_um,
            wavelengths_m=spec_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "spc-spec_long",
                "final_sampling_lam0" => 0.1,
            ),
            description="Compact SPC spec-long model over 15% band",
        ),
        "full_spc_spec_long" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=spec_wavelengths_um,
            wavelengths_m=spec_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "spc-spec_long",
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
            ),
            description="Full SPC spec-long model over 15% band without error maps",
        ),
        "compact_spc_wide" => (
            func=wfirst_phaseb_compact,
            output_dim=128,
            wavelengths_um=wide_wavelengths_um,
            wavelengths_m=wide_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "spc-wide",
                "final_sampling_lam0" => 0.1,
            ),
            description="Compact SPC wide-field model over 10% band",
        ),
        "full_spc_wide" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=wide_wavelengths_um,
            wavelengths_m=wide_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "spc-wide",
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
            ),
            description="Full SPC wide-field model over 10% band without error maps",
        ),
        "compact_hlc_source_offset" => (
            func=wfirst_phaseb_compact,
            output_dim=128,
            wavelengths_um=hlc_wavelengths_um,
            wavelengths_m=hlc_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "hlc",
                "use_hlc_dm_patterns" => 1,
                "final_sampling_lam0" => 0.1,
                "source_x_offset" => 3.0,
            ),
            description="Compact HLC model with nonzero lambda-D source offset",
        ),
        "full_hlc_no_field_stop" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=hlc_wavelengths_um,
            wavelengths_m=hlc_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "hlc",
                "use_hlc_dm_patterns" => 1,
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
                "use_field_stop" => 0,
            ),
            description="Full HLC model without field stop",
        ),
        "full_spc_spec_long_no_pupil_mask" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=spec_wavelengths_um,
            wavelengths_m=spec_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "spc-spec_long",
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
                "use_pupil_mask" => 0,
            ),
            description="Full SPC spec-long model without pupil mask",
        ),
        "full_none" => (
            func=wfirst_phaseb,
            output_dim=128,
            wavelengths_um=hlc_wavelengths_um,
            wavelengths_m=hlc_wavelengths_m,
            passvalue=Dict(
                "cor_type" => "none",
                "final_sampling_lam0" => 0.1,
                "use_errors" => 0,
                "use_fpm" => 0,
            ),
            description="Full pupil-only model without coronagraph elements",
        ),
    )
end

function prepare_phaseb_models(case::NamedTuple; data_root::AbstractString=phaseb_default_data_root())
    passvalue = merge(case.passvalue, Dict("data_dir" => String(data_root)))
    cor_type = String(case.passvalue["cor_type"])
    shared = nothing
    if startswith(cor_type, "hlc")
        shared = prepare_phaseb_shared_assets(data_root, case.wavelengths_m)
    end
    models = map(enumerate(zip(case.wavelengths_um, case.wavelengths_m))) do (i, (λum, λm))
        if shared === nothing
            return prepare_model(
                Symbol("wfirst_" * string(i)),
                case.func,
                λum,
                case.output_dim;
                PASSVALUE=passvalue,
                pool_size=1,
            )
        end
        return prepare_model(
            Symbol("wfirst_" * string(i)),
            case.func,
            λum,
            case.output_dim;
            PASSVALUE=passvalue,
            assets=prepare_asset_pool(() -> PhaseBPreparedAssets(shared, _nearest_occulter(shared, λm), PhaseBModelWorkspace(case.output_dim)); pool_size=1),
            pool_size=1,
        )
    end
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
