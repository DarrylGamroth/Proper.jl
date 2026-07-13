using Test
using Random

function _warmed_gpu_qphase_alloc(wf, z, ctx, sync!)
    prop_qphase(wf, z, ctx)
    sync!()
    prop_qphase(wf, z, ctx)
    sync!()
    return @allocated begin
        prop_qphase(wf, z, ctx)
        sync!()
    end
end

function _warmed_gpu_ptp_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.PLANAR
    prop_ptp(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.PLANAR
    prop_ptp(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        prop_ptp(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_wts_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.PLANAR
    prop_wts(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.PLANAR
    prop_wts(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        prop_wts(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_stw_alloc(wf, dz, ctx, sync!)
    wf.reference_surface = Proper.SPHERICAL
    prop_stw(wf, dz, ctx)
    sync!()
    wf.reference_surface = Proper.SPHERICAL
    prop_stw(wf, dz, ctx)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.SPHERICAL
        prop_stw(wf, dz, ctx)
        sync!()
    end
end

function _warmed_gpu_end_real_alloc(out, wf, sync!)
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf)
    sync!()
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        Proper.prop_end!(out, wf)
        sync!()
    end
end

function _warmed_gpu_end_complex_alloc(out, wf, sync!)
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf; noabs=true)
    sync!()
    wf.reference_surface = Proper.PLANAR
    Proper.prop_end!(out, wf; noabs=true)
    sync!()
    return @allocated begin
        wf.reference_surface = Proper.PLANAR
        Proper.prop_end!(out, wf; noabs=true)
        sync!()
    end
end

function _warmed_gpu_8th_mask_alloc(mask, wf, sync!)
    prop_8th_order_mask!(
        mask,
        wf,
        3f0;
        circular=true,
        min_transmission=0.1f0,
        max_transmission=0.9f0,
    )
    sync!()
    prop_8th_order_mask!(
        mask,
        wf,
        3f0;
        circular=true,
        min_transmission=0.1f0,
        max_transmission=0.9f0,
    )
    sync!()
    return @allocated begin
        prop_8th_order_mask!(
            mask,
            wf,
            3f0;
            circular=true,
            min_transmission=0.1f0,
            max_transmission=0.9f0,
        )
        sync!()
    end
end

const GPU_WARM_QPHASE_ALLOC_MAX = 12_288
const GPU_WARM_PTP_ALLOC_MAX = 12_288
const GPU_WARM_WTS_ALLOC_MAX = 12_288
const GPU_WARM_STW_ALLOC_MAX = 8_192
const GPU_WARM_END_REAL_ALLOC_MAX = 16_384
const GPU_WARM_END_COMPLEX_ALLOC_MAX = 16_384
const GPU_WARM_8TH_MASK_ALLOC_MAX = 16_384

function _gpu_carrier_phase_smoke!(device_array, sync!; atol, rtol)
    T = Float32
    λ = T(500e-9)
    n = 16
    envelope_wf = Proper.WaveFront(
        device_array(fill(complex(one(T)), n, n)),
        λ,
        T(1e-3),
        zero(T),
        one(T),
    )
    carrier_wf = Proper.WaveFront(
        device_array(fill(complex(one(T)), n, n)),
        λ,
        T(1e-3),
        zero(T),
        one(T),
    )
    envelope_ctx = RunContext(
        typeof(envelope_wf.field),
        envelope_wf.workspace;
        carrier_phase=EnvelopeOnly(),
    )
    carrier_ctx = RunContext(
        typeof(carrier_wf.field),
        carrier_wf.workspace;
        carrier_phase=TrackCarrierPhase(),
    )

    prop_ptp(envelope_wf, λ / T(4), envelope_ctx)
    prop_ptp(carrier_wf, λ / T(4), carrier_ctx)
    sync!()
    envelope = Array(envelope_wf.field)
    carrier = Array(carrier_wf.field)
    @test isapprox(carrier, envelope .* Complex{T}(im); atol=atol, rtol=rtol)
    @test isapprox(abs2.(carrier), abs2.(envelope); atol=atol, rtol=rtol)
    return nothing
end

function _gpu_propagate_to_fft_scratch!(wf, ctx, dz, propagation::Symbol)
    if propagation === :ptp
        wf.reference_surface = Proper.PLANAR
        prop_ptp(wf, dz, ctx)
    elseif propagation === :wts
        wf.reference_surface = Proper.PLANAR
        prop_wts(wf, dz, ctx)
    elseif propagation === :stw
        wf.reference_surface = Proper.SPHERICAL
        prop_stw(wf, dz, ctx)
    else
        throw(ArgumentError("unsupported propagation path: $propagation"))
    end
    return wf
end

function _gpu_8th_mask_scratch_alias_regression!(
    device_array,
    sync!,
    ::Type{T},
    propagation::Symbol;
    n::Int=16,
) where {T<:AbstractFloat}
    field = device_array(fill(complex(one(T), T(0.125)), n, n))
    wf_alias = Proper.WaveFront(field, T(500e-9), T(1e-3), zero(T), one(T))
    ctx = RunContext(wf_alias)
    _gpu_propagate_to_fft_scratch!(wf_alias, ctx, T(0.01), propagation)
    sync!()

    # Planned FFT paths intentionally publish the workspace buffer as the
    # current field. The mask reduction must not use that live field as temp.
    @test wf_alias.field === wf_alias.workspace.fft.scratch
    fft_scratch = wf_alias.workspace.fft.scratch
    forward_plan = wf_alias.workspace.fft.forward_plan
    backward_plan = wf_alias.workspace.fft.backward_plan
    snapshot = copy(wf_alias.field)
    sync!()

    wf_oracle = Proper.WaveFront(
        copy(snapshot),
        wf_alias.wavelength_m,
        wf_alias.sampling_m,
        wf_alias.z_m,
        wf_alias.beam_diameter_m,
    )
    wf_oracle.current_fratio = wf_alias.current_fratio
    @test wf_oracle.field !== wf_oracle.workspace.fft.scratch
    @test wf_oracle.workspace.fft.scratch !== wf_alias.workspace.fft.scratch

    mask_alias = similar(wf_alias.field, T, n, n)
    mask_oracle = similar(wf_oracle.field, T, n, n)
    prop_8th_order_mask!(
        mask_alias,
        wf_alias,
        T(3);
        circular=true,
        min_transmission=T(0.1),
        max_transmission=T(0.9),
    )
    prop_8th_order_mask!(
        mask_oracle,
        wf_oracle,
        T(3);
        circular=true,
        min_transmission=T(0.1),
        max_transmission=T(0.9),
    )
    sync!()

    # Both paths start from one device snapshot and run the same mask kernels;
    # a small epsilon allowance covers backend contraction without masking any
    # corruption of the live field by the reduction temporary.
    tolerance = T(16) * eps(T)
    @test isapprox(Array(mask_alias), Array(mask_oracle); atol=tolerance, rtol=tolerance)
    @test isapprox(Array(wf_alias.field), Array(wf_oracle.field); atol=tolerance, rtol=tolerance)
    @test wf_alias.workspace.fft.scratch === fft_scratch
    @test wf_alias.workspace.fft.forward_plan === forward_plan
    @test wf_alias.workspace.fft.backward_plan === backward_plan
    @test wf_alias.workspace.fft.plans_valid
    return nothing
end

function _gpu_map_apply_smoke!(
    wf_gpu,
    gpu_real_matrix,
    expected_backend_array_type,
    sync!,
)
    n = size(wf_gpu.field, 1)
    scale_map_cpu = fill(Float32(1.25), n, n)
    phase_map_cpu = fill(Float32(1f-9), n, n)
    scale_map_gpu = copyto!(similar(gpu_real_matrix, size(scale_map_cpu)...), scale_map_cpu)
    phase_map_gpu = copyto!(similar(gpu_real_matrix, size(phase_map_cpu)...), phase_map_cpu)

    wf_mul_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_mul_gpu.field, ComplexF32(1))
    wf_mul_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_multiply(wf_mul_gpu, scale_map_gpu)
    sync!()
    prop_multiply(wf_mul_cpu, scale_map_cpu)
    @test isapprox(Array(wf_mul_gpu.field), wf_mul_cpu.field; atol=1f-6, rtol=1f-6)

    wf_div_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_div_gpu.field, ComplexF32(2))
    wf_div_cpu = Proper.WaveFront(fill(ComplexF32(2), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_divide(wf_div_gpu, scale_map_gpu)
    sync!()
    prop_divide(wf_div_cpu, scale_map_cpu)
    @test isapprox(Array(wf_div_gpu.field), wf_div_cpu.field; atol=1f-6, rtol=1f-6)

    wf_phase_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_phase_gpu.field, ComplexF32(1))
    wf_phase_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_add_phase(wf_phase_gpu, phase_map_gpu)
    sync!()
    prop_add_phase(wf_phase_cpu, phase_map_cpu)
    @test isapprox(Array(wf_phase_gpu.field), wf_phase_cpu.field; atol=1f-6, rtol=1f-6)

    mktempdir() do dir
        f = joinpath(dir, "map.fits")
        Proper.prop_fits_write(f, scale_map_cpu; HEADER=Dict("PIXSIZE" => wf_gpu.sampling_m))

        map_read = prop_readmap(wf_gpu, f; SAMPLING=wf_gpu.sampling_m)
        sync!()
        @test map_read isa expected_backend_array_type

        wf_err_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
        fill!(wf_err_gpu.field, ComplexF32(1))
        wf_err_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
        prop_errormap(wf_err_gpu, f; SAMPLING=wf_gpu.sampling_m, WAVEFRONT=true)
        sync!()
        prop_errormap(wf_err_cpu, f; SAMPLING=wf_gpu.sampling_m, WAVEFRONT=true)
        @test isapprox(Array(wf_err_gpu.field), wf_err_cpu.field; atol=1f-5, rtol=1f-5)
    end

    wf_psd_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_psd_gpu.field, ComplexF32(1))
    dmap_gpu = prop_psd_errormap(wf_psd_gpu, 1e-18, 10.0, 3.0; no_apply=true, rng=Random.MersenneTwister(7))
    sync!()
    @test dmap_gpu isa expected_backend_array_type

    wf_psd_apply_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_psd_apply_gpu.field, ComplexF32(1))
    wf_psd_apply_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_psd_errormap(wf_psd_apply_gpu, 1e-18, 10.0, 3.0; MIRROR=true, rng=Random.MersenneTwister(11))
    sync!()
    prop_psd_errormap(wf_psd_apply_cpu, 1e-18, 10.0, 3.0; MIRROR=true, rng=Random.MersenneTwister(11))
    @test isapprox(Array(wf_psd_apply_gpu.field), wf_psd_apply_cpu.field; atol=3f-4, rtol=1f-3)
end

function _gpu_centered_wavefront_helpers_smoke!(wf_gpu, expected_backend_array_type, sync!)
    for (ny, nx) in ((7, 7), (8, 8), (7, 8), (8, 7))
        nsamples = ny * nx
        raw = reshape(
            ComplexF32.(1:nsamples) .+ im .* reverse(ComplexF32.(1:nsamples)),
            ny,
            nx,
        )
        field_gpu = copyto!(similar(wf_gpu.field, ComplexF32, ny, nx), raw)
        wf = Proper.WaveFront(
            field_gpu,
            wf_gpu.wavelength_m,
            wf_gpu.sampling_m,
            wf_gpu.z_m,
            wf_gpu.beam_diameter_m,
        )

        centered_field = prop_get_wavefront(wf)
        centered_amplitude = prop_get_amplitude(wf)
        centered_phase = prop_get_phase(wf)
        sync!()

        @test centered_field isa expected_backend_array_type
        @test centered_amplitude isa expected_backend_array_type
        @test centered_phase isa expected_backend_array_type
        @test Array(centered_field) == prop_shift_center(raw)
        @test isapprox(Array(centered_amplitude), prop_shift_center(abs.(raw)); atol=2f-6, rtol=2f-6)
        @test isapprox(Array(centered_phase), prop_shift_center(angle.(raw)); atol=2f-6, rtol=2f-6)

        fill!(centered_field, 0)
        sync!()
        @test Array(wf.field) == raw

        centered_addition = reverse(raw; dims=1)
        addition_gpu = copyto!(similar(wf.field), centered_addition)
        scalar = ComplexF32(2, -3)
        @test prop_add_wavefront(wf, addition_gpu) === wf
        @test prop_add_wavefront(wf, scalar) === wf
        sync!()
        @test Array(wf.field) == raw .+ prop_shift_center(centered_addition; inverse=true) .+ scalar
    end

    ny, nx = 7, 8
    parent_cpu = reshape(
        ComplexF32.(1:(2 * ny * nx)) .+ im .* reverse(ComplexF32.(1:(2 * ny * nx))),
        2 * ny,
        nx,
    )
    view_cpu = parent_cpu[1:2:end, :]
    parent_gpu = copyto!(similar(wf_gpu.field, ComplexF32, size(parent_cpu)...), parent_cpu)
    view_gpu = @view parent_gpu[1:2:end, :]
    shifted_gpu = prop_shift_center(view_gpu)
    @test Proper.backend_style(typeof(view_gpu)) == Proper.backend_style(typeof(parent_gpu))
    @test Proper._shift_center_exec_style(
        typeof(shifted_gpu),
        typeof(view_gpu),
        ny,
        nx,
    ) isa Proper.ShiftKAStyle
    sync!()
    @test shifted_gpu isa expected_backend_array_type
    @test Array(shifted_gpu) == prop_shift_center(view_cpu)

    roundtrip_gpu = prop_shift_center(shifted_gpu; inverse=true)
    sync!()
    @test Array(roundtrip_gpu) == view_cpu

    output_parent_gpu = similar(parent_gpu)
    fill!(output_parent_gpu, 0)
    output_view_gpu = @view output_parent_gpu[1:2:end, :]
    @test Proper._shift_center_exec_style(
        typeof(output_view_gpu),
        typeof(view_gpu),
        ny,
        nx,
    ) isa Proper.ShiftKAStyle
    prop_shift_center!(output_view_gpu, view_gpu)
    sync!()
    output_parent_cpu = Array(output_parent_gpu)
    @test output_parent_cpu[1:2:end, :] == prop_shift_center(view_cpu)
    @test all(iszero, @view output_parent_cpu[2:2:end, :])

    field_cpu = reverse(view_cpu; dims=2)
    field_gpu = copyto!(similar(wf_gpu.field, ComplexF32, ny, nx), field_cpu)
    wf_view_add = Proper.WaveFront(
        field_gpu,
        wf_gpu.wavelength_m,
        wf_gpu.sampling_m,
        wf_gpu.z_m,
        wf_gpu.beam_diameter_m,
    )
    @test prop_add_wavefront(wf_view_add, view_gpu) === wf_view_add
    sync!()
    @test Array(wf_view_add.field) == field_cpu .+ prop_shift_center(view_cpu; inverse=true)
    return nothing
end

function _gpu_direct_dm_smoke!(wf_gpu, gpu_real_matrix, expected_backend_array_type, sync!)
    n = size(wf_gpu.field, 1)
    dm_cpu = fill(Float32(1f-9), n, n)
    dm_gpu = copyto!(similar(gpu_real_matrix, size(dm_cpu)...), dm_cpu)

    wf_dm_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_dm_gpu.field, ComplexF32(1))
    wf_dm_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    @test prop_dm(wf_dm_gpu, dm_gpu) === wf_dm_gpu
    sync!()
    prop_dm(wf_dm_cpu, dm_cpu)
    @test wf_dm_gpu.field isa expected_backend_array_type
    @test isapprox(Array(wf_dm_gpu.field), wf_dm_cpu.field; atol=1f-5, rtol=1f-5)

    wf_mirror_gpu = Proper.WaveFront(similar(wf_gpu.field), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    fill!(wf_mirror_gpu.field, ComplexF32(1))
    wf_mirror_cpu = Proper.WaveFront(fill(ComplexF32(1), n, n), wf_gpu.wavelength_m, wf_gpu.sampling_m, wf_gpu.z_m, wf_gpu.beam_diameter_m)
    prop_dm(wf_mirror_gpu, dm_gpu; mirror=true)
    sync!()
    prop_dm(wf_mirror_cpu, dm_cpu; mirror=true)
    @test wf_mirror_gpu.field isa expected_backend_array_type
    @test isapprox(Array(wf_mirror_gpu.field), wf_mirror_cpu.field; atol=1f-5, rtol=1f-5)
end

@testset "R2 trait-driven routing" begin
    @testset "Style dispatch for cubic convolution" begin
        struct TestInterpStyle <: Proper.InterpStyle end
        Proper._cubic_sample(::TestInterpStyle, a::AbstractMatrix, y::Real, x::Real) = 123.0

        a = reshape(collect(1.0:16.0), 4, 4)
        out = Proper.prop_cubic_conv(TestInterpStyle(), a, [1.0, 2.0], [1.0, 2.0]; grid=true)
        @test size(out) == (2, 2)
        @test all(out .== 123.0)
    end

    @testset "Context-routed interpolation kernels" begin
        wf = prop_begin(1.0, 500e-9, 16)
        dmap = rand(TEST_RNG, Float32, 16, 16)
        opts = Proper.ResampleMapOptions(wf, wf.sampling_m, 8.0, 8.0)
        ctx = RunContext(typeof(dmap))

        out_ctx = similar(dmap, Float64, size(wf.field)...)
        out_def = similar(dmap, Float64, size(wf.field)...)
        Proper.prop_resamplemap!(out_ctx, wf, dmap, opts, ctx)
        Proper.prop_resamplemap!(out_def, wf, dmap, opts)
        @test isapprox(out_ctx, out_def; atol=0, rtol=0)

        img = rand(TEST_RNG, Float32, 16, 16)
        r_ctx = prop_rotate(img, 12.0, ctx)
        r_def = prop_rotate(img, 12.0)
        @test isapprox(r_ctx, r_def; atol=0, rtol=0)

        m_ctx = prop_magnify(img, 1.2, 16, ctx; QUICK=true)
        m_def = prop_magnify(img, 1.2, 16; QUICK=true)
        @test isapprox(m_ctx, m_def; atol=1e-6, rtol=1e-6)
    end

    @testset "Explicit backend boundary for hot output paths" begin
        wf = prop_begin(1.0, 500e-9, 16)
        cpu_out = zeros(Float64, 16, 16)
        @test prop_end!(cpu_out, wf) === cpu_out
    end

    @testset "KA interpolation pilot parity on large CPU arrays" begin
        n = 256
        expected_ka = Threads.nthreads() > 1 && n * n >= Proper.cpu_interp_ka_min_elems()
        @test Proper.cpu_interp_ka_min_elems(1) == 16_384
        @test Proper.cpu_interp_ka_min_elems(16) == 32_768
        @test Proper.ka_cubic_grid_enabled(Matrix{Float32}, n, n) == expected_ka
        @test Proper.ka_rotate_enabled(Matrix{Float32}, n, n) == expected_ka
        Proper.with_cpu_inner_kernel_parallelism(false) do
            @test !Proper.ka_cubic_grid_enabled(Matrix{Float32}, n, n)
            @test !Proper.ka_rotate_enabled(Matrix{Float32}, n, n)
        end

        a = reshape(collect(Float32, 1:(n * n)), n, n)
        x = collect(Float32, 1:n)
        y = collect(Float32, 1:n)
        out_loop = similar(a)
        out_ka = similar(a)
        Proper._prop_cubic_conv_grid_loop!(out_loop, Proper.CubicInterpStyle(), a, x, y)
        Proper.ka_cubic_conv_grid!(out_ka, a, x, y)
        @test isapprox(out_ka, out_loop; atol=0, rtol=0)

        opts_cubic = Proper.RotateOptions(a, pairs((; CUBIC=true)))
        rc_loop = similar(a)
        rc_ka = similar(a)
        c = cos(-deg2rad(9.0))
        s = sin(-deg2rad(9.0))
        Proper._prop_rotate_cubic!(Proper.CubicInterpStyle(), rc_loop, a, c, s, opts_cubic)
        Proper.ka_rotate_cubic!(rc_ka, a, c, s, opts_cubic.cx, opts_cubic.cy, opts_cubic.sx, opts_cubic.sy)
        @test isapprox(rc_ka, rc_loop; atol=0, rtol=0)

        opts_linear = Proper.RotateOptions(a, pairs((; METH="linear")))
        rl_loop = similar(a)
        rl_ka = similar(a)
        Proper._prop_rotate_linear!(rl_loop, a, c, s, opts_linear)
        Proper.ka_rotate_linear!(rl_ka, a, c, s, opts_linear.cx, opts_linear.cy, opts_linear.sx, opts_linear.sy)
        @test isapprox(rl_ka, rl_loop; atol=0, rtol=0)
    end

    @testset "KA geometry and sampling parity on large CPU arrays" begin
        n = 256
        @test !Proper.ka_geometry_enabled(Matrix{Float32}, n, n)
        @test !Proper.ka_sampling_enabled(Matrix{Float32}, n, n)

        wf = prop_begin(1.0, 500e-9, n)
        RT = real(eltype(wf.field))

        rect_loop = zeros(RT, n, n)
        rect_ka = similar(rect_loop)
        Proper.prop_rectangle!(rect_loop, wf, 0.4, 0.2, 0.03, -0.05; ROTATION=22.0, NORM=true)
        dx = RT(prop_get_sampling(wf))
        beamrad = RT(prop_get_beamradius(wf))
        pr = beamrad / dx
        θ = RT(deg2rad(22.0))
        cθ = cos(θ)
        sθ = sin(θ)
        xcp = RT(n ÷ 2) + RT(0.03) * pr
        ycp = RT(n ÷ 2) - RT(0.05) * pr
        xrp = RT(0.4) * pr / RT(2)
        yrp = RT(0.2) * pr / RT(2)
        xp = (-xrp, -xrp, xrp, xrp)
        yp = (-yrp, yrp, yrp, -yrp)
        xbox = ntuple(i -> xp[i] * cθ - yp[i] * sθ + xcp, 4)
        ybox = ntuple(i -> xp[i] * sθ + yp[i] * cθ + ycp, 4)
        minx = max(0, floor(Int, minimum(xbox) - one(RT)))
        maxx = min(n - 1, ceil(Int, maximum(xbox) + one(RT)))
        miny = max(0, floor(Int, minimum(ybox) - one(RT)))
        maxy = min(n - 1, ceil(Int, maximum(ybox) + one(RT)))
        Proper.ka_rectangle_mask!(rect_ka, xcp, ycp, xrp, yrp, cθ, sθ, minx, maxx, miny, maxy; nsub=Proper.antialias_subsampling())
        @test isapprox(rect_ka, rect_loop; atol=1e-12, rtol=1e-12)

        ellipse_loop = zeros(RT, n, n)
        ellipse_ka = similar(ellipse_loop)
        Proper.prop_ellipse!(ellipse_loop, wf, 0.35, 0.25, 0.05, -0.03; ROTATION=13.0, DARK=true, NORM=true)
        beamrad_pix = RT(prop_get_beamradius(wf)) / dx
        xcenter_pix = RT(n ÷ 2) + RT(0.05) * beamrad_pix
        ycenter_pix = RT(n ÷ 2) - RT(0.03) * beamrad_pix
        xrad_pix = RT(0.35) * beamrad_pix
        yrad_pix = RT(0.25) * beamrad_pix
        t = RT(deg2rad(13.0))
        sint = sin(t)
        cost = cos(t)
        delx = inv(xrad_pix)
        dely = inv(yrad_pix)
        drx = delx * cost - dely * sint
        dry = delx * sint + dely * cost
        dr = max(abs(drx), abs(dry))
        Proper.ka_ellipse_mask!(ellipse_ka, xcenter_pix, ycenter_pix, xrad_pix, yrad_pix, sint, cost, one(RT) + dr, one(RT) - dr, RT(1 + 1e-10); dark=true, nsub=Proper.antialias_subsampling())
        @test isapprox(ellipse_ka, ellipse_loop; atol=1e-12, rtol=1e-12)

        xverts = RT[-0.20, 0.12, 0.28, -0.08]
        yverts = RT[-0.18, -0.22, 0.19, 0.25]
        ipoly_loop = zeros(RT, n, n)
        ipoly_ka = similar(ipoly_loop)
        Proper.prop_irregular_polygon!(ipoly_loop, wf, xverts, yverts; NORM=true)
        xv = copy(xverts) .* beamrad
        yv = copy(yverts) .* beamrad
        Proper.ka_irregular_polygon_mask!(ipoly_ka, xv, yv, n ÷ 2, n ÷ 2, dx; nsub=Proper.antialias_subsampling())
        @test isapprox(ipoly_ka, ipoly_loop; atol=1e-12, rtol=1e-12)

        round_loop = Proper.prop_rounded_rectangle(wf, 0.05, 0.3, 0.2, 0.01, -0.02)
        round_ka = similar(round_loop)
        Proper.ka_rounded_rectangle_mask!(round_ka, dx, RT(0.01), RT(-0.02), RT(0.05), RT(0.3) / RT(2), RT(0.2) / RT(2))
        @test isapprox(round_ka, round_loop; atol=1e-12, rtol=1e-12)

        img = rand(TEST_RNG, Float32, n, n)
        pix_loop = Proper._prop_pixellate_factor(img, 2)
        pix_ka = similar(pix_loop)
        Proper.ka_pixellate!(pix_ka, img, 2)
        @test isapprox(pix_ka, pix_loop; atol=0, rtol=0)

        mag = Float32(1.35)
        n_out = 128
        szoom_loop = prop_szoom(img, mag, n_out)
        szoom_ka = similar(szoom_loop)
        table_loop = Matrix{Float32}(undef, n_out, Proper.SZOOM_K)
        table_ka = similar(table_loop)
        Proper._fill_szoom_table_loop!(table_loop, mag)
        Proper.ka_szoom_table!(table_ka, mag, n_out, Proper.SZOOM_K, Proper.SZOOM_DK)
        @test isapprox(table_ka, table_loop; atol=0, rtol=0)
        Proper.ka_szoom_apply!(szoom_ka, img, table_ka, mag)
        @test isapprox(szoom_ka, szoom_loop; atol=1e-6, rtol=1e-6)

        rect_img = rand(TEST_RNG, Float32, n, n + 16)
        rect_loop = prop_szoom(rect_img, mag; NOX=120, NOY=96)
        rect_ka = similar(rect_loop)
        tablex_loop = Matrix{Float32}(undef, 120, Proper.SZOOM_K)
        tabley_loop = Matrix{Float32}(undef, 96, Proper.SZOOM_K)
        tablex_ka = similar(tablex_loop)
        tabley_ka = similar(tabley_loop)
        Proper._fill_szoom_table_loop!(tablex_loop, mag)
        Proper._fill_szoom_table_loop!(tabley_loop, mag)
        Proper.ka_szoom_table!(tablex_ka, mag, 120, Proper.SZOOM_K, Proper.SZOOM_DK)
        Proper.ka_szoom_table!(tabley_ka, mag, 96, Proper.SZOOM_K, Proper.SZOOM_DK)
        @test isapprox(tablex_ka, tablex_loop; atol=0, rtol=0)
        @test isapprox(tabley_ka, tabley_loop; atol=0, rtol=0)
        Proper.ka_szoom_apply!(rect_ka, rect_img, tablex_ka, tabley_ka, mag)
        @test isapprox(rect_ka, rect_loop; atol=1e-6, rtol=1e-6)
    end

    @testset "Context-routed propagation kernels" begin
        wf1 = prop_begin(1.0, 500e-9, 32)
        wf2 = prop_begin(1.0, 500e-9, 32)
        ctx = RunContext(typeof(wf2.field))

        prop_propagate(wf1, 0.05)
        prop_propagate(wf2, 0.05, ctx)

        @test wf1.z_m == wf2.z_m
        @test isapprox(wf1.field, wf2.field; atol=0, rtol=0)

        lens1 = prop_begin(1.0, 500e-9, 32)
        lens2 = prop_begin(1.0, 500e-9, 32)
        lensctx = RunContext(lens2)
        prop_lens(lens1, 0.35)
        prop_lens(lens2, 0.35, lensctx)
        @test lens1.reference_surface === lens2.reference_surface
        @test lens1.beam_type_old === lens2.beam_type_old
        @test lens1.propagator_type === lens2.propagator_type
        @test isapprox(lens1.field, lens2.field; atol=0, rtol=0)
    end

    @testset "Optional CUDA smoke (no scalar indexing)" begin
        cuda_ready = false
        try
            @eval using CUDA
            cuda_ready = CUDA.functional()
        catch
            cuda_ready = false
        end

        if cuda_ready
            CUDA.allowscalar(false)
            @test Base.get_extension(Proper, :ProperCUDAExt) !== nothing

            a = CUDA.rand(Float32, 16, 16)
            ctx = RunContext(typeof(a))
            @test ctx.backend isa Proper.CUDABackend
            @test ctx.fft isa Proper.CUFFTStyle
            @test ctx.interp isa Proper.CubicInterpStyle
            @test Proper.ka_cubic_grid_enabled(typeof(a), 16, 16)
            @test Proper.ka_rotate_enabled(typeof(a), 16, 16)
            @test Proper.ka_geometry_enabled(typeof(a), 16, 16)
            @test Proper.ka_sampling_enabled(typeof(a), 16, 16)
            @test Proper.qphase_workspace(ctx).xphase isa CUDA.CuArray
            @test Proper.qphase_workspace(ctx).yphase isa CUDA.CuArray
            @test !Proper.ka_separable_phase_enabled(typeof(a), 16, 16)
            @test Proper.ka_separable_phase_enabled(typeof(a), 256, 256)
            @test !Proper.ka_separable_qphase_enabled(typeof(a), 256, 256)
            @test Proper.ka_separable_qphase_enabled(CUDA.CuMatrix{ComplexF64}, 128, 128)
            _gpu_carrier_phase_smoke!(CUDA.CuArray, CUDA.synchronize; atol=1f-5, rtol=1f-5)

            m = prop_magnify(a, 1.1, 16, ctx; QUICK=true)
            r = prop_rotate(a, 5.0, ctx)
            s = prop_szoom(a, 1.1, 16)
            p = Proper._prop_pixellate_factor(a, 2)
            wf_resample = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            ropts = Proper.ResampleMapOptions(wf_resample, wf_resample.sampling_m, 8f0, 8f0)
            res = similar(a)
            Proper.prop_resamplemap!(res, wf_resample, a, ropts, ctx)
            @test_throws ArgumentError Proper.prop_cubic_conv(a, 1.5f0, 1.5f0)
            img_cpu = reshape(Float32.(1:256), 16, 16)
            xcoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), 1, :), 8, 1)
            ycoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), :, 1), 1, 8)
            coord_ref = Proper.prop_cubic_conv(img_cpu, xcoords_cpu, ycoords_cpu; grid=false)
            coord_gpu = Proper.prop_cubic_conv(CUDA.CuArray(img_cpu), xcoords_cpu, ycoords_cpu; grid=false)
            @test coord_gpu isa CUDA.CuArray
            @test isapprox(Array(coord_gpu), coord_ref; atol=1f-5, rtol=1f-5)
            @test size(m) == (16, 16)
            @test size(r) == (16, 16)
            @test size(s) == (16, 16)
            @test size(p) == (8, 8)
            @test size(res) == (16, 16)
            @test Proper.interp_workspace(ctx).xcoords isa CUDA.CuArray
            @test Proper.interp_workspace(ctx).ycoords isa CUDA.CuArray

            wf = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf, 0.25f0, ctx)
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            fws = Proper.fft_workspace(ctx)
            @test fws.scratch isa CUDA.CuArray
            rho2 = Proper.ensure_rho2_map!(fws, 16, 16, 1f-3)
            @test rho2 isa CUDA.CuArray
            @test fws.forward_plan !== nothing
            @test fws.backward_plan !== nothing
            pfft = fws.forward_plan
            pbfft = fws.backward_plan
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            @test fws.forward_plan === pfft
            @test fws.backward_plan === pbfft
            @test Proper.ensure_rho2_map!(fws, 16, 16, 1f-3) === rho2
            prop_circular_aperture(wf, 2.5f-4)
            @test wf.workspace.mask.mask isa CUDA.CuArray
            wf_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf_ref, 0.25f0)
            wf_ref.reference_surface = Proper.PLANAR
            prop_ptp(wf_ref, 0.01f0)
            prop_circular_aperture(wf_ref, 2.5f-4)
            @test isapprox(Array(wf.field), wf_ref.field; atol=1f-5, rtol=1f-5)
            wf_mask = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            wf_mask_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            mask = prop_8th_order_mask(wf_mask, 3f0; circular=true)
            mask_ref = prop_8th_order_mask(wf_mask_ref, 3f0; circular=true)
            @test mask isa CUDA.CuArray
            @test isapprox(Array(mask), mask_ref; atol=1f-5, rtol=1f-5)
            @test isapprox(Array(wf_mask.field), wf_mask_ref.field; atol=1f-5, rtol=1f-5)
            @test wf_mask.workspace.mask.reduction_scratch isa CUDA.CuVector{ComplexF32}
            rect = prop_rectangle(wf, 5f-4, 4f-4)
            round = prop_rounded_rectangle(wf, 2f-4, 5f-4, 4f-4)
            out, sampling = prop_end(wf)
            @test_throws ArgumentError Proper.prop_end!(zeros(Float32, 16, 16), wf)
            out_ref, sampling_ref = prop_end(wf_ref)
            @test size(rect) == (16, 16)
            @test size(round) == (16, 16)
            @test size(out) == (16, 16)
            @test sampling == wf.sampling_m
            @test sampling == sampling_ref
            @test isapprox(Array(out), out_ref; atol=1f-5, rtol=1f-5)
            wf_alloc = Proper.WaveFront(CUDA.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            out_alloc = similar(wf_alloc.field, Float32, 16, 16)
            @test _warmed_gpu_qphase_alloc(wf_alloc, 0.25f0, ctx, CUDA.synchronize) <= GPU_WARM_QPHASE_ALLOC_MAX
            @test _warmed_gpu_ptp_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_PTP_ALLOC_MAX
            @test _warmed_gpu_wts_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_WTS_ALLOC_MAX
            @test _warmed_gpu_stw_alloc(wf_alloc, 0.01f0, ctx, CUDA.synchronize) <= GPU_WARM_STW_ALLOC_MAX
            @test _warmed_gpu_end_real_alloc(out_alloc, wf_alloc, CUDA.synchronize) <= GPU_WARM_END_REAL_ALLOC_MAX
            @test _warmed_gpu_end_complex_alloc(similar(wf_alloc.field), wf_alloc, CUDA.synchronize) <= GPU_WARM_END_COMPLEX_ALLOC_MAX
            _gpu_direct_dm_smoke!(wf_alloc, CUDA.zeros(Float32, 16, 16), CUDA.CuArray, CUDA.synchronize)
            _gpu_map_apply_smoke!(wf_alloc, CUDA.zeros(Float32, 16, 16), CUDA.CuArray, CUDA.synchronize)
            _gpu_centered_wavefront_helpers_smoke!(wf_alloc, CUDA.CuArray, CUDA.synchronize)
            @testset "8th-order mask after planned FFT" begin
                for T in (Float32, Float64), propagation in (:ptp, :wts, :stw)
                    @testset "$propagation $T" begin
                        n = T === Float32 && propagation === :ptp ? 256 : 16
                        _gpu_8th_mask_scratch_alias_regression!(
                            CUDA.CuArray,
                            CUDA.synchronize,
                            T,
                            propagation;
                            n,
                        )
                    end
                end
            end
        else
            @test true
        end
    end

    @testset "Optional AMDGPU smoke (no scalar indexing)" begin
        amdgpu_ready = false
        try
            @eval using AMDGPU
            amdgpu_ready = AMDGPU.functional() && AMDGPU.functional(:rocfft)
        catch
            amdgpu_ready = false
        end

        if amdgpu_ready
            AMDGPU.allowscalar(false)
            @test Base.get_extension(Proper, :ProperAMDGPUExt) !== nothing

            a = AMDGPU.rand(Float32, 16, 16)
            ctx = RunContext(typeof(a))
            @test ctx.backend isa Proper.AMDGPUBackend
            @test ctx.fft isa Proper.ROCFFTStyle
            @test ctx.interp isa Proper.CubicInterpStyle
            @test Proper.ka_cubic_grid_enabled(typeof(a), 16, 16)
            @test Proper.ka_geometry_enabled(typeof(a), 16, 16)
            @test Proper.ka_sampling_enabled(typeof(a), 16, 16)
            @test Proper.qphase_workspace(ctx).xphase isa AMDGPU.ROCArray
            @test Proper.qphase_workspace(ctx).yphase isa AMDGPU.ROCArray
            @test !Proper.ka_separable_phase_enabled(typeof(a), 16, 16)
            @test Proper.ka_separable_phase_enabled(typeof(a), 256, 256)
            @test !Proper.ka_separable_qphase_enabled(typeof(a), 256, 256)
            @test Proper.ka_separable_qphase_enabled(AMDGPU.ROCMatrix{ComplexF64}, 128, 128)
            _gpu_carrier_phase_smoke!(AMDGPU.ROCArray, AMDGPU.synchronize; atol=3f-4, rtol=1f-3)

            promoted_c = 10.0
            promoted_wf_gpu = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-6, 0f0, 1f0)
            promoted_ctx = RunContext(promoted_wf_gpu)
            promoted_k = pi / (promoted_wf_gpu.wavelength_m * promoted_c)
            promoted_ref = Matrix{ComplexF32}(undef, 16, 16)
            for j in axes(promoted_ref, 2), i in axes(promoted_ref, 1)
                x = Float32(Proper._ka_shifted_index_0based(j - 1, 16)) * promoted_wf_gpu.sampling_m
                y = Float32(Proper._ka_shifted_index_0based(i - 1, 16)) * promoted_wf_gpu.sampling_m
                promoted_ref[i, j] = cis(promoted_k * (x * x + y * y))
            end
            prop_qphase(promoted_wf_gpu, promoted_c, promoted_ctx)
            AMDGPU.synchronize()
            @test isapprox(Array(promoted_wf_gpu.field), promoted_ref; atol=3f-6, rtol=3f-6)
            @test isempty(Proper.qphase_workspace(promoted_ctx).xphase)
            @test isempty(Proper.qphase_workspace(promoted_ctx).yphase)

            p = Proper._prop_pixellate_factor(a, 2)
            @test_throws ArgumentError Proper.prop_cubic_conv(a, 1.5f0, 1.5f0)
            img_cpu = reshape(Float32.(1:256), 16, 16)
            xcoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), 1, :), 8, 1)
            ycoords_cpu = repeat(reshape(collect(Float32, range(3, 13; length=8)), :, 1), 1, 8)
            coord_ref = Proper.prop_cubic_conv(img_cpu, xcoords_cpu, ycoords_cpu; grid=false)
            coord_gpu = Proper.prop_cubic_conv(AMDGPU.ROCArray(img_cpu), xcoords_cpu, ycoords_cpu; grid=false)
            @test coord_gpu isa AMDGPU.ROCArray
            @test isapprox(Array(coord_gpu), coord_ref; atol=3f-4, rtol=1f-3)
            @test size(p) == (8, 8)

            wf = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf, 0.25f0, ctx)
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            fws = Proper.fft_workspace(ctx)
            @test fws.scratch isa AMDGPU.ROCArray
            rho2 = Proper.ensure_rho2_map!(fws, 16, 16, 1f-3)
            @test rho2 isa AMDGPU.ROCArray
            @test fws.forward_plan !== nothing
            @test fws.backward_plan !== nothing
            pfft = fws.forward_plan
            pbfft = fws.backward_plan
            wf.reference_surface = Proper.PLANAR
            prop_ptp(wf, 0.01f0, ctx)
            @test fws.forward_plan === pfft
            @test fws.backward_plan === pbfft
            @test Proper.ensure_rho2_map!(fws, 16, 16, 1f-3) === rho2
            wf_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            prop_qphase(wf_ref, 0.25f0)
            wf_ref.reference_surface = Proper.PLANAR
            prop_ptp(wf_ref, 0.01f0)
            @test isapprox(Array(wf.field), wf_ref.field; atol=3f-4, rtol=1f-3)

            phase_n = 256
            phase_rng = MersenneTwister(20260712)
            phase_input = rand(phase_rng, ComplexF64, phase_n, phase_n)
            wf_phase_ref = Proper.WaveFront(copy(phase_input), 500e-9, 1e-6, 0.0, 1.0)
            wf_phase_gpu = Proper.WaveFront(AMDGPU.ROCArray(phase_input), 500e-9, 1e-6, 0.0, 1.0)
            ctx_phase = RunContext(wf_phase_gpu)
            prop_qphase(wf_phase_ref, 10.0)
            prop_qphase(wf_phase_gpu, 10.0, ctx_phase)
            AMDGPU.synchronize()
            @test isapprox(Array(wf_phase_gpu.field), wf_phase_ref.field; atol=3e-10, rtol=3e-10)
            @test length(Proper.qphase_workspace(ctx_phase).xphase) == phase_n
            @test length(Proper.qphase_workspace(ctx_phase).yphase) == phase_n

            ptp_input = rand(phase_rng, ComplexF32, phase_n, phase_n)
            wf_ptp_ref = Proper.WaveFront(copy(ptp_input), 500f-9, 1f-4, 0f0, 1f0)
            wf_ptp_gpu = Proper.WaveFront(AMDGPU.ROCArray(ptp_input), 500f-9, 1f-4, 0f0, 1f0)
            ctx_ptp = RunContext(wf_ptp_gpu)
            wf_ptp_ref.reference_surface = Proper.PLANAR
            wf_ptp_gpu.reference_surface = Proper.PLANAR
            prop_ptp(wf_ptp_ref, 0.01f0)
            prop_ptp(wf_ptp_gpu, 0.01f0, ctx_ptp)
            AMDGPU.synchronize()
            @test isapprox(Array(wf_ptp_gpu.field), wf_ptp_ref.field; atol=3f-4, rtol=1f-3)
            @test length(Proper.qphase_workspace(ctx_ptp).xphase) == phase_n
            @test length(Proper.qphase_workspace(ctx_ptp).yphase) == phase_n

            wf_mask = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            wf_mask_ref = Proper.WaveFront(fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            mask = prop_8th_order_mask(wf_mask, 3f0; circular=true)
            mask_ref = prop_8th_order_mask(wf_mask_ref, 3f0; circular=true)
            @test mask isa AMDGPU.ROCArray
            @test isapprox(Array(mask), mask_ref; atol=3f-4, rtol=1f-3)
            @test isapprox(Array(wf_mask.field), wf_mask_ref.field; atol=3f-4, rtol=1f-3)

            mask_rng = MersenneTwister(20260713)
            mask_cases = (
                (;
                    dims=(15, 17),
                    circular=false,
                    elliptical=nothing,
                    y_axis=true,
                    min_transmission=0.05f0,
                    max_transmission=0.8f0,
                ),
                (;
                    dims=(16, 18),
                    circular=true,
                    elliptical=nothing,
                    y_axis=false,
                    min_transmission=0.1f0,
                    max_transmission=0.9f0,
                ),
                (;
                    dims=(15, 17),
                    circular=false,
                    elliptical=1.4f0,
                    y_axis=false,
                    min_transmission=0.0f0,
                    max_transmission=0.7f0,
                ),
            )
            for mask_case in mask_cases
                mask_input = rand(mask_rng, ComplexF32, mask_case.dims...)
                wf_mask_cpu = Proper.WaveFront(copy(mask_input), 500f-9, 1f-3, 0f0, 1f0)
                wf_mask_gpu = Proper.WaveFront(AMDGPU.ROCArray(mask_input), 500f-9, 1f-3, 0f0, 1f0)
                mask_cpu = zeros(Float32, mask_case.dims...)
                mask_gpu = AMDGPU.zeros(Float32, mask_case.dims...)

                prop_8th_order_mask!(
                    mask_cpu,
                    wf_mask_cpu,
                    3f0;
                    circular=mask_case.circular,
                    elliptical=mask_case.elliptical,
                    y_axis=mask_case.y_axis,
                    min_transmission=mask_case.min_transmission,
                    max_transmission=mask_case.max_transmission,
                )
                prop_8th_order_mask!(
                    mask_gpu,
                    wf_mask_gpu,
                    3f0;
                    circular=mask_case.circular,
                    elliptical=mask_case.elliptical,
                    y_axis=mask_case.y_axis,
                    min_transmission=mask_case.min_transmission,
                    max_transmission=mask_case.max_transmission,
                )
                AMDGPU.synchronize()

                @test isapprox(Array(mask_gpu), mask_cpu; atol=3f-5, rtol=3f-5)
                @test isapprox(Array(wf_mask_gpu.field), wf_mask_cpu.field; atol=3f-5, rtol=3f-5)
                @test wf_mask_gpu.workspace.mask.reduction_scratch isa AMDGPU.ROCArray
            end

            wf_mask_alloc = Proper.WaveFront(
                AMDGPU.fill(ComplexF32(1), 32, 34),
                500f-9,
                1f-3,
                0f0,
                1f0,
            )
            mask_alloc = AMDGPU.zeros(Float32, 32, 34)
            mask_alloc_bytes = _warmed_gpu_8th_mask_alloc(
                mask_alloc,
                wf_mask_alloc,
                AMDGPU.synchronize,
            )
            @test mask_alloc_bytes <= GPU_WARM_8TH_MASK_ALLOC_MAX
            reduction_scratch = wf_mask_alloc.workspace.mask.reduction_scratch
            @test reduction_scratch isa AMDGPU.ROCVector{ComplexF32}
            @test length(reduction_scratch) == Proper._mask_extrema_temp_length(length(mask_alloc))
            prop_8th_order_mask!(mask_alloc, wf_mask_alloc, 3f0; circular=true)
            AMDGPU.synchronize()
            @test wf_mask_alloc.workspace.mask.reduction_scratch === reduction_scratch

            prop_circular_aperture(wf, 2.5f-4)
            @test wf.workspace.mask.mask isa AMDGPU.ROCArray
            prop_circular_aperture(wf_ref, 2.5f-4)
            out, sampling = prop_end(wf)
            @test_throws ArgumentError Proper.prop_end!(zeros(Float32, 16, 16), wf)
            out_ref, sampling_ref = prop_end(wf_ref)
            @test size(out) == (16, 16)
            @test sampling == wf.sampling_m
            @test sampling == sampling_ref
            @test Array(out) isa Matrix{Float32}
            wf_alloc = Proper.WaveFront(AMDGPU.fill(ComplexF32(1), 16, 16), 500f-9, 1f-3, 0f0, 1f0)
            out_alloc = similar(wf_alloc.field, Float32, 16, 16)
            @test _warmed_gpu_qphase_alloc(wf_alloc, 0.25f0, ctx, AMDGPU.synchronize) <= GPU_WARM_QPHASE_ALLOC_MAX
            @test _warmed_gpu_ptp_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_PTP_ALLOC_MAX
            @test _warmed_gpu_wts_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_WTS_ALLOC_MAX
            @test _warmed_gpu_stw_alloc(wf_alloc, 0.01f0, ctx, AMDGPU.synchronize) <= GPU_WARM_STW_ALLOC_MAX
            @test _warmed_gpu_end_real_alloc(out_alloc, wf_alloc, AMDGPU.synchronize) <= GPU_WARM_END_REAL_ALLOC_MAX
            @test _warmed_gpu_end_complex_alloc(similar(wf_alloc.field), wf_alloc, AMDGPU.synchronize) <= GPU_WARM_END_COMPLEX_ALLOC_MAX
            _gpu_direct_dm_smoke!(wf_alloc, AMDGPU.zeros(Float32, 16, 16), AMDGPU.ROCArray, AMDGPU.synchronize)
            _gpu_centered_wavefront_helpers_smoke!(wf_alloc, AMDGPU.ROCArray, AMDGPU.synchronize)
            @testset "8th-order mask after planned FFT" begin
                for T in (Float32, Float64), propagation in (:ptp, :wts, :stw)
                    @testset "$propagation $T" begin
                        n = T === Float32 && propagation === :ptp ? 256 : 16
                        _gpu_8th_mask_scratch_alias_regression!(
                            AMDGPU.ROCArray,
                            AMDGPU.synchronize,
                            T,
                            propagation;
                            n,
                        )
                    end
                end
            end
        else
            @test true
        end
    end
end
