using Test
using Random

@testset "CPU KernelAbstractions helper kernels" begin
    rng = MersenneTwister(20260428)

    field = fill(ComplexF32(1 + 0im), 8, 8)
    mask = zeros(Float32, 8, 8)
    mask[4:6, 4:6] .= 1
    shifted_mask = copy(field)
    shifted_mask_ref = copy(field)
    @test Proper.ka_apply_shifted_mask!(shifted_mask, mask) === shifted_mask
    Proper._apply_shifted_mask_loop!(shifted_mask_ref, mask)
    @test shifted_mask == shifted_mask_ref

    shifted_mask_invert = copy(field)
    shifted_mask_invert_ref = copy(field)
    @test Proper.ka_apply_shifted_mask!(shifted_mask_invert, mask; invert=true) === shifted_mask_invert
    Proper._apply_shifted_mask_loop!(shifted_mask_invert_ref, mask; invert=true)
    @test shifted_mask_invert == shifted_mask_invert_ref

    @testset "Fused 8th-order mask normalization and shifted application" begin
        cases = (((5, 7), false), ((6, 8), false), ((3, 5), true), ((1, 1), true))
        for T in (Float32, Float64), (dims, degenerate) in cases
            raw = degenerate ? fill(T(0.37), dims) : rand(rng, T, dims...)
            min_transmission = T(0.1)
            max_transmission = T(0.9)

            squared_mask = raw .* raw
            reference_mask = copy(squared_mask)
            squared_min = minimum(reference_mask)
            squared_max = maximum(reference_mask)
            reference_mask .-= squared_min
            squared_range = maximum(reference_mask)
            if squared_range > 0
                reference_mask ./= squared_range
            end
            reference_mask .*= max_transmission - min_transmission
            reference_mask .+= min_transmission
            reference_mask .= sqrt.(reference_mask)

            candidate_mask = copy(raw)
            candidate_field = fill(complex(one(T), T(0.25)), dims)
            reference_field = candidate_field .* prop_shift_center(reference_mask)

            @test Proper.ka_normalize_apply_8th_mask!(
                candidate_mask,
                candidate_field,
                squared_min,
                squared_max,
                min_transmission,
                max_transmission,
            ) === candidate_mask
            @test candidate_mask == reference_mask
            @test candidate_field == reference_field
        end

        @test_throws ArgumentError Proper.ka_normalize_apply_8th_mask!(
            zeros(Float32, 3, 5),
            zeros(ComplexF32, 5, 3),
            0.0f0,
            1.0f0,
            0.0f0,
            1.0f0,
        )
    end

    @test Proper.ka_apply_centered_circle!(copy(field), 9.0f0, 4.0f0, 6.25f0; nsub=2) isa Matrix{ComplexF32}
    @test Proper.ka_apply_centered_circle!(copy(field), 9.0f0, 4.0f0, 6.25f0; dark=true, invert=true, nsub=2) isa Matrix{ComplexF32}
    @test Proper.ka_apply_shifted_circle!(copy(field), 0.25f0, -0.5f0, 2.5f0, 9.0f0, 4.0f0, 6.25f0; nsub=2) isa Matrix{ComplexF32}
    @test Proper.ka_apply_shifted_ellipse!(
        copy(field),
        4.0f0,
        4.0f0,
        2.5f0,
        1.5f0,
        sin(0.25f0),
        cos(0.25f0),
        1.1f0,
        0.9f0,
        1.0f0;
        minx_pix=1,
        maxx_pix=6,
        miny_pix=1,
        maxy_pix=6,
        nsub=2,
    ) isa Matrix{ComplexF32}

    src = reshape(ComplexF32.(1:64), 8, 8)
    outc = zeros(ComplexF32, 4, 4)
    @test Proper.ka_copy_shifted_complex!(outc, src, 2, 2) === outc
    outi = zeros(Float32, 4, 4)
    @test Proper.ka_copy_shifted_intensity!(outi, src, 2, 2) === outi

    shifted = zeros(Float32, 8, 8)
    @test Proper.ka_shift_copy!(shifted, Float32.(reshape(1:64, 8, 8)), 2, 3) === shifted
    @test shifted == circshift(Float32.(reshape(1:64, 8, 8)), (2, 3))

    qfield = fill(ComplexF32(1 + 0im), 8, 8)
    @test Proper.ka_apply_qphase!(qfield, 0.5f0, 0.1f0) === qfield
    qfield_ref = Matrix{ComplexF32}(undef, 8, 8)
    for j in axes(qfield_ref, 2), i in axes(qfield_ref, 1)
        x0 = Float32(Proper._ka_shifted_index_0based(j - 1, 8)) * 0.1f0
        y0 = Float32(Proper._ka_shifted_index_0based(i - 1, 8)) * 0.1f0
        qfield_ref[i, j] = cis(0.5f0 * (x0 * x0 + y0 * y0))
    end
    @test qfield ≈ qfield_ref
    @test Proper.ka_apply_frequency_phase!(qfield, 0.25f0, 0.1f0) === qfield
    @test !Proper.ka_separable_phase_enabled(Matrix{ComplexF32}, 256, 256)
    @test !Proper.ka_separable_qphase_enabled(Matrix{ComplexF64}, 256, 256)

    promoted_field = fill(ComplexF32(1), 5, 5)
    promoted_wf = Proper.WaveFront(promoted_field, 500f-9, 1f-6, 0f0, 1f0)
    promoted_ctx = RunContext(promoted_wf)
    promoted_c = 10.0
    promoted_k = pi / (promoted_wf.wavelength_m * promoted_c)
    @test promoted_k isa Float64
    promoted_ref = Matrix{ComplexF32}(undef, size(promoted_field))
    for j in axes(promoted_ref, 2), i in axes(promoted_ref, 1)
        x = Float32(Proper._ka_shifted_index_0based(j - 1, size(promoted_ref, 2))) * promoted_wf.sampling_m
        y = Float32(Proper._ka_shifted_index_0based(i - 1, size(promoted_ref, 1))) * promoted_wf.sampling_m
        promoted_ref[i, j] = cis(promoted_k * (x * x + y * y))
    end
    @test Proper._prop_qphase_ka!(promoted_wf, promoted_c, promoted_ctx) === promoted_wf
    @test promoted_wf.field ≈ promoted_ref

    sep_ws = Proper.QPhaseWorkspace(Float32)
    sep_k = 0.5f0
    sep_sy = 0.2f0
    sep_sx = 0.1f0
    for sep_scale in (1.0f0, inv(7.0f0))
        sep_field = fill(ComplexF32(1 + 0.25im), 7, 9)
        sep_ref = copy(sep_field)
        xphase_ref = Vector{ComplexF32}(undef, size(sep_ref, 2))
        yphase_ref = Vector{ComplexF32}(undef, size(sep_ref, 1))
        Proper._fill_fft_order_axis_phase!(xphase_ref, length(xphase_ref), sep_sx, sep_k)
        Proper._fill_fft_order_axis_phase!(yphase_ref, length(yphase_ref), sep_sy, sep_k)
        Proper._apply_separable_phase!(sep_ref, xphase_ref, yphase_ref, sep_scale)
        @test Proper._apply_separable_quadratic_phase_ka!(
            sep_field,
            sep_k,
            sep_sy,
            sep_sx,
            sep_ws,
            sep_scale,
        ) === sep_field
        @test sep_field ≈ sep_ref
        @test length(sep_ws.xphase) == size(sep_field, 2)
        @test length(sep_ws.yphase) == size(sep_field, 1)
    end

    xaxis = zeros(Float32, 5)
    yaxis = zeros(Float32, 3)
    @test Proper.ka_fill_affine_axis!(xaxis, 2.0f0, 0.5f0, 1.0f0) === xaxis
    @test xaxis ≈ Float32[0.0, 0.5, 1.0, 1.5, 2.0]
    @test Proper.ka_fill_affine_axes!(xaxis, yaxis, 2.0f0, 1.0f0, 0.5f0, 2.0f0, 1.0f0, -1.0f0) === (xaxis, yaxis)

    rho2 = zeros(Float32, 8, 8)
    @test Proper.ka_fill_fft_order_rho2!(rho2, 0.25f0) === rho2
    @test rho2[1, 1] == 0
    @test Proper.ka_scale_field!(qfield, 0.5f0) === qfield

    image = rand(rng, Float32, 8, 8)
    xval = collect(Float32, range(2, 7; length=4))
    yval = collect(Float32, range(2, 7; length=4))
    cubic = zeros(Float32, 4, 4)
    @test Proper.ka_cubic_conv_grid!(cubic, image, xval, yval) === cubic
    cubic_ref = similar(cubic)
    Proper._prop_cubic_conv_grid_loop!(cubic_ref, Proper.CubicInterpStyle(), image, xval, yval)
    @test cubic ≈ cubic_ref

    xgrid = repeat(reshape(xval, 1, :), length(yval), 1)
    ygrid = repeat(reshape(yval, :, 1), 1, length(xval))
    cubic_coord = zeros(Float32, size(xgrid))
    @test Proper.ka_cubic_conv_coordinate_grid!(cubic_coord, image, xgrid, ygrid) === cubic_coord
    cubic_coord_ref = Proper._prop_cubic_conv(Proper.CubicInterpStyle(), Proper.PointwiseTopology(), image, xgrid, ygrid)
    @test cubic_coord ≈ cubic_coord_ref

    cubic_coord_mut = similar(cubic_coord)
    @test Proper.prop_cubic_conv_coordinate_grid!(cubic_coord_mut, Proper.CubicInterpStyle(), image, xgrid, ygrid) === cubic_coord_mut
    @test cubic_coord_mut ≈ cubic_coord_ref

    rotlin = zeros(Float32, 8, 8)
    rotlin_ref = similar(rotlin)
    lin_opts = Proper.RotateOptions(image, pairs((; METH="linear")))
    c = cos(0.1f0)
    s = sin(0.1f0)
    @test Proper.ka_rotate_linear!(rotlin, image, c, s, lin_opts.cx, lin_opts.cy, lin_opts.sx, lin_opts.sy) === rotlin
    Proper._prop_rotate_linear!(rotlin_ref, image, c, s, lin_opts)
    @test rotlin ≈ rotlin_ref

    rotcub = zeros(Float32, 8, 8)
    rotcub_ref = similar(rotcub)
    cubic_opts = Proper.RotateOptions(image, pairs((; CUBIC=true)))
    @test Proper.ka_rotate_cubic!(rotcub, image, c, s, cubic_opts.cx, cubic_opts.cy, cubic_opts.sx, cubic_opts.sy) === rotcub
    Proper._prop_rotate_cubic!(Proper.CubicInterpStyle(), rotcub_ref, image, c, s, cubic_opts)
    @test rotcub ≈ rotcub_ref

    rect = zeros(Float32, 8, 8)
    @test Proper.ka_rectangle_mask!(rect, 4.0f0, 4.0f0, 2.0f0, 1.5f0, cos(0.2f0), sin(0.2f0), 1, 6, 1, 6; nsub=2) === rect
    @test 0 < sum(rect) < length(rect)

    ell = zeros(Float32, 8, 8)
    @test Proper.ka_ellipse_mask!(ell, 4.0f0, 4.0f0, 2.5f0, 1.5f0, sin(0.2f0), cos(0.2f0), 1.1f0, 0.9f0, 1.0f0; nsub=2) === ell
    ell_ref = Proper.prop_ellipse(prop_begin(1.0, 500e-9, 8), 2.5f-3, 1.5f-3, 0.0f0, 0.0f0; ROTATION=rad2deg(0.2f0))
    @test 0 < sum(ell) < length(ell)
    @test size(ell_ref) == size(ell)

    poly = zeros(Float32, 8, 8)
    xv = Float32[-0.2, 0.2, 0.2, -0.2]
    yv = Float32[-0.2, -0.2, 0.2, 0.2]
    @test Proper.ka_irregular_polygon_mask!(poly, xv, yv, 4, 4, 0.1f0; nsub=2) === poly
    @test 0 < sum(poly) < length(poly)

    rounded = zeros(Float32, 8, 8)
    @test Proper.ka_rounded_rectangle_mask!(rounded, 0.1f0, 0.0f0, 0.0f0, 0.1f0, 0.25f0, 0.2f0) === rounded
    @test 0 < sum(rounded) < length(rounded)

    table = zeros(Float32, 8, Proper.SZOOM_K)
    @test Proper.ka_szoom_table!(table, 1.2f0, 8, Proper.SZOOM_K, Proper.SZOOM_DK) === table
    tablex = zeros(Float32, 8, Proper.SZOOM_K)
    tabley = zeros(Float32, 6, Proper.SZOOM_K)
    @test Proper.ka_szoom_tables!(tablex, tabley, 1.2f0, Proper.SZOOM_K, Proper.SZOOM_DK) === (tablex, tabley)
    zout = zeros(Float32, 8, 8)
    @test Proper.ka_szoom_apply!(zout, image, table, 1.2f0) === zout
    @test zout ≈ prop_szoom(image, 1.2f0, 8)

    zrect = zeros(Float32, 6, 8)
    @test Proper.ka_szoom_apply!(zrect, image, tablex, tabley, 1.2f0) === zrect
    @test zrect ≈ prop_szoom(image, 1.2f0; NOX=8, NOY=6)

    pix = zeros(Float32, 4, 4)
    @test Proper.ka_pixellate!(pix, image, 2) === pix
    @test pix ≈ Proper._prop_pixellate_factor(image, 2)
end
