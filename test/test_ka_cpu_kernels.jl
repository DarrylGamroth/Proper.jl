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
