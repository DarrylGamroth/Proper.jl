using Test

function exercise_gpu_coronagraph_diagnostics!(mod, device_array, synchronize)
    input = reshape(ComplexF32.(1:256), 16, 16)
    wf = Proper.WaveFront(device_array(input), 500f-9, 1f-3, 0f0, 1f0)
    diagnostics = mod.CoronagraphDiagnostics(wf)

    @test Proper.same_backend_style(
        typeof(diagnostics.after_occulter_amplitude),
        typeof(wf.field),
    )
    @test Proper.same_backend_style(
        typeof(diagnostics.before_lyot_amplitude),
        typeof(wf.field),
    )
    @test eltype(diagnostics.after_occulter_amplitude) === Float32
    @test eltype(diagnostics.before_lyot_amplitude) === Float32
    mod._capture_after_occulter!(diagnostics, wf)
    mod._capture_before_lyot!(diagnostics, wf)
    synchronize()

    expected = prop_shift_center(abs.(input))
    @test Array(diagnostics.after_occulter_amplitude) == expected
    @test Array(diagnostics.before_lyot_amplitude) == expected
    return nothing
end

@testset "Example ports smoke" begin
    exdir = joinpath(@__DIR__, "..", "examples")
    examples = (
        ("simple_prescription.jl", :simple_prescription),
        ("simple_telescope.jl", :simple_telescope),
        ("hubble_simple.jl", :hubble_simple),
        ("microscope.jl", :microscope),
        ("run_coronagraph.jl", :run_coronagraph),
        ("run_coronagraph_dm.jl", :run_coronagraph_dm),
        ("run_occulter.jl", :run_occulter),
        ("talbot.jl", :talbot),
        ("talbot_correct.jl", :talbot_correct),
        ("psdtest.jl", :psdtest),
        ("multi_example.jl", :multi_example),
        ("migration_dm_fits.jl", :migration_dm_fits_demo),
        ("wfirst_phaseb_reference.jl", :wfirst_phaseb_reference_demo),
    )

    for (file, sym) in examples
        path = joinpath(exdir, file)
        mod = load_example_module(path)
        @test isdefined(mod, sym)
        @test !occursin(r"(?m)^using Plots$", read(path, String))
    end
end

@testset "Coronagraph diagnostics preserve final numerics" begin
    mod = load_example_module(joinpath(@__DIR__, "..", "examples", "coronagraph.jl"))
    n = 32
    diam = 0.1
    f_lens = 24 * diam

    function initial_wavefront()
        wf = prop_begin(diam, 0.55e-6, n; beam_diam_fraction=0.3)
        prop_circular_aperture(wf, diam / 2)
        prop_define_entrance(wf)
        return wf
    end

    reference_wf = initial_wavefront()
    diagnostic_wf = initial_wavefront()
    diagnostics = mod.CoronagraphDiagnostics(diagnostic_wf)

    @test mod.coronagraph(reference_wf, f_lens, :eighth_order, diam) === reference_wf
    @test mod.coronagraph(
        diagnostic_wf,
        f_lens,
        :eighth_order,
        diam;
        diagnostics,
    ) === diagnostic_wf

    reference_psf, reference_sampling = prop_end(reference_wf)
    diagnostic_psf, diagnostic_sampling = prop_end(diagnostic_wf)
    @test diagnostic_psf == reference_psf
    @test diagnostic_sampling == reference_sampling

    for amplitude in (
        diagnostics.after_occulter_amplitude,
        diagnostics.before_lyot_amplitude,
    )
        @test size(amplitude) == (n, n)
        @test all(isfinite, amplitude)
        @test all(>=(0), amplitude)
        @test maximum(amplitude) > 0
    end

    bad_diagnostics = mod.CoronagraphDiagnostics(
        zeros(n - 1, n - 1),
        zeros(n - 1, n - 1),
    )
    @test_throws ArgumentError mod.coronagraph(
        initial_wavefront(),
        f_lens,
        :eighth_order,
        diam;
        diagnostics=bad_diagnostics,
    )

    cuda_ready = false
    try
        @eval using CUDA
        cuda_ready = CUDA.functional()
    catch
        cuda_ready = false
    end
    if cuda_ready
        CUDA.allowscalar(false)
        exercise_gpu_coronagraph_diagnostics!(mod, CUDA.CuArray, CUDA.synchronize)
    else
        @test true
    end

    amdgpu_ready = false
    try
        @eval using AMDGPU
        amdgpu_ready = AMDGPU.functional()
    catch
        amdgpu_ready = false
    end
    if amdgpu_ready
        AMDGPU.allowscalar(false)
        exercise_gpu_coronagraph_diagnostics!(mod, AMDGPU.ROCArray, AMDGPU.synchronize)
    else
        @test true
    end
end

@testset "Migration DM/FITS example executes" begin
    mod = load_example_module(joinpath(@__DIR__, "..", "examples", "migration_dm_fits.jl"))
    psf, sampling = mod.migration_dm_fits_demo(; gridsize=32)
    @test size(psf) == (32, 32)
    @test sampling > 0
    @test maximum(psf) > 0
end

@testset "Multi example executes with upstream actuator-space DM map" begin
    mod = load_example_module(joinpath(@__DIR__, "..", "examples", "multi_example.jl"))
    field, sampling = mod.multi_example(0.55e-6, 32, Dict("use_dm" => true, "dm" => zeros(48, 48)))
    @test size(field) == (32, 32)
    @test sampling > 0
end

@testset "Multi-run example regressions execute" begin
    testmulti1_mod = load_example_module(joinpath(@__DIR__, "..", "examples", "testmulti1.jl"))
    psf = testmulti1_mod.testmulti1(; gridsize=32, npsf=32, nlambda=2)
    @test size(psf) == (32, 32)
    @test all(isfinite, psf)

    testmulti2_mod = load_example_module(joinpath(@__DIR__, "..", "examples", "testmulti2.jl"))
    fields, samplings = testmulti2_mod.testmulti2(; gridsize=32, npatterns=2)
    @test size(fields) == (32, 32, 2)
    @test all(isfinite, fields)
    @test all(>(0), samplings)
end

@testset "Example native keyword APIs remain compatible with PASSVALUE" begin
    exdir = joinpath(@__DIR__, "..", "examples")

    run_occulter_mod = load_example_module(joinpath(exdir, "run_occulter.jl"))
    native_solid, native_sampling = run_occulter_mod.run_occulter(0.55e-6, 32; occulter=:solid)
    compat_solid, compat_sampling = run_occulter_mod.run_occulter(0.55e-6, 32, Dict("occulter_type" => "SOLID"))
    @test native_solid ≈ compat_solid
    @test native_sampling == compat_sampling

    talbot_mod = load_example_module(joinpath(exdir, "talbot.jl"))
    native_field, native_talbot_sampling = talbot_mod.talbot(0.5e-6, 32; diam=0.1, period=0.04, dist=0.0)
    compat_field, compat_talbot_sampling = talbot_mod.talbot(0.5e-6, 32, Dict("diam" => 0.1, "period" => 0.04, "dist" => 0.0))
    @test native_field ≈ compat_field
    @test native_talbot_sampling == compat_talbot_sampling
end
