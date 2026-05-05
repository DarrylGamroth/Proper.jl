using Test
using Random
using Proper

function _contract_wavefront_prescription(λm, n; scale=1.0)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.4 * scale)
    prop_define_entrance(wf)
    return wf
end

function _contract_tuple_prescription(λm, n; scale=1.0)
    wf = prop_begin(1.0, λm, n)
    prop_rectangle(wf, 0.2 * scale, 0.3 * scale)
    return prop_end(wf)
end

function _contract_passvalue_prescription(λm, n, passvalue; scale=1.0)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.2 + 0.05 * float(passvalue) * scale)
    prop_define_entrance(wf)
    return prop_end(wf)
end

function _contract_keyword_passvalue_prescription(λm, n; radius=0.4, scale=1.0)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, radius * scale)
    prop_define_entrance(wf)
    return prop_end(wf)
end

@testset "API contract smoke tests" begin
    wf = prop_begin(1.0, 550e-9, 32)
    @test wf isa WaveFront

    out_intensity, sampling = prop_end(wf)
    @test out_intensity isa AbstractMatrix
    @test size(out_intensity) == (32, 32)
    @test sampling isa Float64

    wf2 = prop_begin(1.0, 550e-9, 32)
    prop_lens(wf2, 5.0)
    prop_propagate(wf2, 0.1)
    @test wf2 isa WaveFront

    dm = zeros(32, 32)
    @test prop_dm(wf2, dm, 16.0, 16.0, 0.05; NO_APPLY=true) isa AbstractMatrix

    img = rand(TEST_RNG, 32, 32)
    @test size(prop_rotate(img, 0.0)) == size(img)
    @test size(prop_magnify(img, 0.5)) == (16, 16)
    @test size(prop_resamplemap(wf2, rand(TEST_RNG, 16, 16), wf2.sampling_m, 8.0, 8.0)) == size(wf2.field)
    @test size(prop_pixellate(img, 0.5, 1.0, 16)) == (16, 16)
    @test size(prop_8th_order_mask(wf2, 3.0)) == size(wf2.field)

    mktempdir() do d
        map = rand(TEST_RNG, 32, 32) .* 1e-9
        mapfile = joinpath(d, "contract_map.fits")
        prop_writemap(map, mapfile; SAMPLING=wf2.sampling_m)
        @test prop_errormap(wf2, mapfile) === wf2

        rng = MersenneTwister(1234)
        dmap = prop_psd_errormap(wf2, 1e-18, 10.0, 3.0; NO_APPLY=true, RNG=rng)
        @test dmap isa AbstractMatrix
        @test size(dmap) == size(wf2.field)
    end
end

@testset "API contract execution shapes" begin
    psf, sampling = prop_run(_contract_wavefront_prescription, 0.55, 16)
    @test psf isa AbstractMatrix
    @test size(psf) == (16, 16)
    @test sampling isa Float64

    psf2, sampling2 = prop_run(_contract_tuple_prescription, 0.55, 16)
    @test psf2 isa AbstractMatrix
    @test size(psf2) == (16, 16)
    @test sampling2 isa Float64

    prepared = prepare_prescription(_contract_passvalue_prescription, 0.55, 16; PASSVALUE=0.0)
    batch = prepare_prescription_batch(prepared; pool_size=2)
    assets = prepare_asset_pool(; pool_size=2) do slot
        return (scale=1.0 + 0.1 * slot,)
    end
    model = prepare_model(prepared; name=:contract_model, assets=assets, pool_size=2)

    ppsf, psampling = prop_run(prepared; PASSVALUE=1.0)
    bpsf, bsampling = prop_run(batch; PASSVALUE=2.0, slot=1)
    mpsf, msampling = prop_run(model; PASSVALUE=3.0, slot=2)
    @test size(ppsf) == (16, 16)
    @test size(bpsf) == (16, 16)
    @test size(mpsf) == (16, 16)
    @test psampling isa Float64
    @test bsampling isa Float64
    @test msampling isa Float64

    stack, samplings = prop_run_multi(_contract_passvalue_prescription, 0.55, 16; PASSVALUE=[0.0, 1.0, 2.0])
    @test stack isa Array{<:Number,3}
    @test size(stack) == (16, 16, 3)
    @test samplings isa Vector{Float64}
    @test length(samplings) == 3

    kpsf, ksampling = prop_run(
        _contract_keyword_passvalue_prescription,
        0.55,
        16;
        PASSVALUE=Dict("RADIUS" => 0.25),
        scale=1.0,
    )
    @test size(kpsf) == (16, 16)
    @test ksampling isa Float64
    @test_throws ArgumentError prop_run(
        _contract_keyword_passvalue_prescription,
        0.55,
        16;
        PASSVALUE=Dict(:radius => 0.25),
        radius=0.30,
    )
end

@testset "API contract keyword compatibility" begin
    img = reshape(collect(1.0:25.0), 5, 5)
    @test prop_rotate(img, 0.0; METH="linear", XSHIFT=1.0, MISSING=-1.0) ==
          prop_rotate(img, 0.0; meth="linear", xshift=1.0, missing=-1.0)

    @test prop_magnify(img, 1.0; QUICK=true) ==
          prop_magnify(img, 1.0; quick=true)

    wf = prop_begin(1.0, 550e-9, 16)
    map = rand(TEST_RNG, 8, 8)
    @test prop_resamplemap(wf, map, wf.sampling_m, 4.0, 4.0, 0.25, -0.5) ==
          prop_resamplemap(wf, map, wf.sampling_m, 4.0, 4.0, 0.25, -0.5)
end

@testset "API contract export surface" begin
    @test Base.isexported(Proper, :prop_run)
    @test Base.isexported(Proper, :prop_run_multi)
    @test Base.isexported(Proper, :prepare_prescription)
    @test Base.isexported(Proper, :prepare_prescription_batch)
    @test Base.isexported(Proper, :prepare_asset_pool)
    @test Base.isexported(Proper, :prepare_model)
    @test Base.isexported(Proper, :prepare_hot_call)
    @test Base.isexported(Proper, :prepared_context)
    @test Base.isexported(Proper, :prop_run_hot)

    @test Base.isexported(Proper, :prop_dftidefs)

    @test !Base.isexported(Proper, :libcconv)
    @test !Base.isexported(Proper, :libcconvthread)
    @test !Base.isexported(Proper, :switch_set)
    @test !Base.isexported(Proper, :DftiErrorMessage)
    @test !isdefined(Proper, :prop_execute_multi)
    @test !isdefined(Proper, :prop_table)
    @test !isdefined(Proper, :prop_fftw)
    @test !isdefined(Proper, :prop_ffti)
    @test !isdefined(Proper, :prop_use_fftw)
    @test !isdefined(Proper, :prop_use_ffti)
    @test !Base.isexported(Proper, :WFIRSTPhaseBProper)
    @test !isdefined(Proper, :WFIRSTPhaseBProper)
end
