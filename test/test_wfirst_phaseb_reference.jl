using Test
using FFTW
using Statistics
using Proper
using Proper.WFIRSTPhaseBProper

@testset "WFIRST Phase B reference helpers" begin
    old_root = data_dir()
    try
        set_data_dir!(joinpath(pwd(), ".cache", "wfirst", "dummy"))
        @test phaseb_default_data_root() == data_dir()
    finally
        set_data_dir!(old_root)
    end

    a = reshape(collect(1.0:16.0), 4, 4)
    @test trim(a, 2) == [6.0 10.0; 7.0 11.0]
    @test size(trim(a, 6)) == (6, 6)

    c = ComplexF64.(reshape(1:16, 4, 4))
    @test ffts(copy(c), -1) ≈ circshift(fft(circshift(c, (-2, -2))) ./ length(c), (2, 2))

    m = mft2(c, 0.1, 2.0, 4, -1)
    @test size(m) == (4, 4)
    let field_in = c, dout = 0.1, D = 2.0, nout = 4, direction = -1
        nfield_in = size(field_in, 2)
        x = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2))
        y = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2))
        u = (collect(0:(nout - 1)) .- (nout ÷ 2)) .* (dout / D)
        v = (collect(0:(nout - 1)) .- (nout ÷ 2)) .* (dout / D)
        xu = x * transpose(u)
        yv = y * transpose(v)
        expxu = (dout / D) .* exp.(-direction * 2π * im .* xu)
        expyv = transpose(exp.(-direction * 2π * im .* yv))
        @test m ≈ expyv * field_in * expxu
    end

    cases = phaseb_case_definitions()
    @test Set(keys(cases)) == Set([
        "compact_hlc",
        "full_hlc",
        "compact_spc_spec_long",
        "full_spc_spec_long",
        "compact_spc_wide",
        "full_spc_wide",
        "compact_hlc_source_offset",
        "full_hlc_no_field_stop",
        "full_spc_spec_long_no_pupil_mask",
        "full_none",
    ])
    @test cases["compact_hlc"].func === wfirst_phaseb_compact
    @test cases["full_hlc"].func === wfirst_phaseb
    @test cases["compact_spc_spec_long"].passvalue["cor_type"] == "spc-spec_long"
    @test cases["full_spc_wide"].passvalue["cor_type"] == "spc-wide"
    @test cases["compact_hlc_source_offset"].passvalue["source_x_offset"] == 3.0
    @test cases["full_hlc_no_field_stop"].passvalue["use_field_stop"] == 0
    @test cases["full_spc_spec_long_no_pupil_mask"].passvalue["use_pupil_mask"] == 0
    @test cases["full_none"].passvalue["cor_type"] == "none"

    sx, sy = Proper.WFIRSTPhaseBProper._source_offset_lambda_over_d((source_x_offset_mas=10.0, source_y_offset=1.5), 0.575e-6, 2.363)
    expected_mas_per_lamd = 0.575e-6 * 360.0 * 3600.0 / (2π * 2.363) * 1000
    @test isapprox(sx, 10.0 / expected_mas_per_lamd; rtol=1e-12)
    @test sy == 1.5

    cfg_hlc = Proper.WFIRSTPhaseBProper._phaseb_config("hlc", 0.575e-6, "/tmp"; compact=false, use_fpm=1)
    @test cfg_hlc.branch == :hlc
    @test cfg_hlc.n_default == 1024
    @test cfg_hlc.n_to_fpm == 2048
    @test occursin("run461_pupil_rotated.fits", cfg_hlc.pupil_file)

    cfg_spc = Proper.WFIRSTPhaseBProper._phaseb_config("spc-ifs_short", 0.66e-6, "/tmp"; compact=false, use_fpm=1)
    @test cfg_spc.branch == :spc
    @test cfg_spc.lambda0_m == 0.66e-6
    @test cfg_spc.n_default == 2048
    @test cfg_spc.n_mft == 1400
    @test occursin("SPM_SPC-20190130.fits", cfg_spc.pupil_mask_file)

    cfg_spc_compact = Proper.WFIRSTPhaseBProper._phaseb_config("spc-wide", 0.825e-6, "/tmp"; compact=true, use_fpm=1)
    @test cfg_spc_compact.branch == :spc
    @test cfg_spc_compact.n_big == 1400
    @test occursin("rotated", cfg_spc_compact.pupil_mask_file)

    cfg_none = Proper.WFIRSTPhaseBProper._phaseb_config("none", 0.575e-6, "/tmp"; compact=false, use_fpm=0)
    @test cfg_none.branch == :none
    @test cfg_none.lyot_stop_file === nothing

    ws = Proper.WFIRSTPhaseBProper.PhaseBModelWorkspace(8)
    @test size(Proper.WFIRSTPhaseBProper.phaseb_field(ws, 1400)) == (1400, 1400)
    @test Proper.WFIRSTPhaseBProper.phaseb_field(ws, 1400) === Proper.WFIRSTPhaseBProper.phaseb_field(ws, 1400)
    @test Proper.WFIRSTPhaseBProper.phaseb_fft_cache(ws, 1400) === Proper.WFIRSTPhaseBProper.phaseb_fft_cache(ws, 1400)

    wf_tilt = prop_begin(1.0, 550e-9, 8)
    field_before_tilt = copy(wf_tilt.field)
    Proper.WFIRSTPhaseBProper._apply_source_offset!(wf_tilt, 6.0, 0.575e-6, 550e-9, 0.1, -0.2)
    @test wf_tilt.field != field_before_tilt

    mktempdir() do d
        old_root = data_dir()
        try
            set_data_dir(d)
            @test data_dir() == abspath(d)

            copied = joinpath(d, "copied")
            copy_here(copied)
            @test isfile(joinpath(copied, "wfirst_phaseb.jl"))
            @test isfile(joinpath(copied, "wfirst_phaseb_compact.jl"))

            copied_examples = joinpath(d, "examples")
            copy_examples_here(copied_examples)
            @test isfile(joinpath(copied_examples, "wfirst_phaseb_reference.jl"))

            polroot = joinpath(d, "coeffs")
            zamp = zeros(Float64, 2, 2, 6, 22)
            zpha = zeros(Float64, 2, 2, 6, 22)
            zamp[1, 1, :, 1] .= 2.0
            prop_fits_write(polroot * "_amp.fits", zamp)
            prop_fits_write(polroot * "_pha.fits", zpha)

            amp, pha = polab(polroot, 550e-9, 6.0, -1)
            @test size(amp) == size(pha)
            @test all(iszero, pha)
            @test isapprox(mean(amp), 2.0; atol=1e-12)

            wf = prop_begin(1.0, 550e-9, 8)
            field_before = copy(wf.field)
            polmap(wf, polroot, 6.0, -1)
            @test wf.field ≈ field_before .* 2
        finally
            set_data_dir!(old_root)
        end
    end
end
