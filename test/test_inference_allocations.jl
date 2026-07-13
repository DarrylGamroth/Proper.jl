using Test
using Random

const JULIA_110_TINY_ALLOC_BYTES = VERSION < v"1.11.0-DEV" ? 64 : 0

@testset "Inference and Allocation Gates" begin
    wf = prop_begin(1.0, 500e-9, 32)
    rng = MersenneTwister(123)
    dmap = rand(rng, Float32, 32, 32)

    opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 16.0, 16.0)
    out = similar(dmap, Float64, size(wf.field)...)
    @test size(Proper.prop_resamplemap!(out, wf, dmap, opts)) == size(wf.field)

    # Warm-up + allocation gate for mutating resample path.
    Proper.prop_resamplemap!(out, wf, dmap, opts)
    resample_alloc = @allocated Proper.prop_resamplemap!(out, wf, dmap, opts)
    @test resample_alloc < 4_000_000

    # Scalar cubic interpolation should remain allocation-free.
    img = rand(rng, Float64, 16, 16)
    _ = prop_cubic_conv(img, 7.3, 8.1)
    @test (@allocated prop_cubic_conv(img, 7.3, 8.1)) <= JULIA_110_TINY_ALLOC_BYTES

    phase_map = rand(rng, Float64, 32, 32)
    scale_map = rand(rng, Float64, 32, 32) .+ 0.5

    wf_phase = prop_begin(1.0, 500e-9, 32)
    prop_add_phase(wf_phase, phase_map)
    add_phase_alloc = @allocated prop_add_phase(wf_phase, phase_map)
    @test add_phase_alloc < 10_000

    wf_mult = prop_begin(1.0, 500e-9, 32)
    prop_multiply(wf_mult, scale_map)
    multiply_alloc = @allocated prop_multiply(wf_mult, scale_map)
    @test multiply_alloc < 10_000

    wf_div = prop_begin(1.0, 500e-9, 32)
    prop_divide(wf_div, scale_map)
    divide_alloc = @allocated prop_divide(wf_div, scale_map)
    @test divide_alloc < 10_000

    mktempdir() do d
        f = joinpath(d, "map.fits")
        prop_fits_write(f, dmap)

        r = @inferred AbstractArray prop_fits_read(f)
        @test size(r) == size(dmap)

        rh, hdr = @inferred Tuple{AbstractArray,Proper.FITSHeader} Proper.prop_fits_read_with_header(f)
        @test size(rh) == size(dmap)
        @test haskey(hdr, "SIMPLE")

    rm = prop_readmap(wf, f; SAMPLING=wf.sampling_m)
    @test size(rm) == size(wf.field)
end

    @test collect(Proper.coordinate_axis(5, 2.0)) == [-4.0, -2.0, 0.0, 2.0, 4.0]
    @test collect(Proper.spatial_frequency_axis(4, 0.5)) == [-1.0, -0.5, 0.0, 0.5]
    @test (@allocated Proper.coordinate_axis(128, 1.0)) <= JULIA_110_TINY_ALLOC_BYTES
    @test (@allocated Proper.spatial_frequency_axis(128, 1.0)) <= JULIA_110_TINY_ALLOC_BYTES

    dummy(λm, n; kwargs...) = prop_begin(1.0, λm, n)
    stack, samplings = prop_run_multi(dummy, 0.55, 16; PASSVALUE=[1, 2, 3])
    @test size(stack) == (16, 16, 3)
    @test eltype(samplings) <: AbstractFloat

    runner = (item, pass, slot) -> (fill(Float32(item + pass), 2, 3), Float32(slot))
    inferred_stack, inferred_samplings = @inferred Proper._prop_run_multi_items(
        Proper.MultiRunThreadedExecStyle(),
        collect(1:3),
        collect(10:12),
        runner,
    )
    @test inferred_stack isa Array{Float32,3}
    @test inferred_samplings isa Vector{Float32}

    asset_type = NamedTuple{(:slot,),Tuple{Int}}
    asset_pool = @inferred prepare_asset_pool((slot, _) -> (slot=slot,), asset_type; pool_size=2)
    @test asset_pool isa PreparedAssetPool
    @test eltype(asset_pool.cache) == Union{Nothing,asset_type}
    asset_model = prepare_model(dummy, 0.55, 16; assets=asset_pool, pool_size=2)
    @test (@inferred Proper.prepared_assets(asset_model, 1)) == (slot=1,)

    @test_throws ArgumentError prepare_asset_pool((slot, _) -> (slot=slot,); pool_size=2)
    @test_throws ArgumentError prepare_asset_pool((slot, _) -> (slot=slot,), Any; pool_size=2)

    bad_asset_pool = prepare_asset_pool((slot, _) -> (slot=string(slot),), asset_type; pool_size=1)
    bad_asset_model = prepare_model(dummy, 0.55, 16; assets=bad_asset_pool, pool_size=1)
    @test_throws ArgumentError Proper.prepared_assets(bad_asset_model, 1)
end
