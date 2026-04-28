using Test
using Random

@testset "CPU KernelAbstractions helper kernels" begin
    rng = MersenneTwister(20260428)

    field = fill(ComplexF32(1 + 0im), 8, 8)
    mask = zeros(Float32, 8, 8)
    mask[4:6, 4:6] .= 1
    @test Proper.ka_apply_shifted_mask!(copy(field), mask) isa Matrix{ComplexF32}
    @test Proper.ka_apply_shifted_mask!(copy(field), mask; invert=true) isa Matrix{ComplexF32}

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
    @test sum(shifted) == sum(Float32.(1:64))

    qfield = fill(ComplexF32(1 + 0im), 8, 8)
    @test Proper.ka_apply_qphase!(qfield, 0.5f0, 0.1f0) === qfield
    @test !all(isone, qfield)
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
    @test all(isfinite, cubic)

    rotlin = zeros(Float32, 8, 8)
    @test Proper.ka_rotate_linear!(rotlin, image, cos(0.1f0), sin(0.1f0), 4.5f0, 4.5f0, 0.0f0, 0.0f0) === rotlin
    rotcub = zeros(Float32, 8, 8)
    @test Proper.ka_rotate_cubic!(rotcub, image, cos(0.1f0), sin(0.1f0), 4.5f0, 4.5f0, 0.0f0, 0.0f0) === rotcub

    rect = zeros(Float32, 8, 8)
    @test Proper.ka_rectangle_mask!(rect, 4.0f0, 4.0f0, 2.0f0, 1.5f0, cos(0.2f0), sin(0.2f0), 1, 6, 1, 6; nsub=2) === rect
    @test 0 < sum(rect) < length(rect)
    ell = zeros(Float32, 8, 8)
    @test Proper.ka_ellipse_mask!(ell, 4.0f0, 4.0f0, 2.5f0, 1.5f0, sin(0.2f0), cos(0.2f0), 1.1f0, 0.9f0, 1.0f0; nsub=2) === ell
    poly = zeros(Float32, 8, 8)
    xv = Float32[-0.2, 0.2, 0.2, -0.2]
    yv = Float32[-0.2, -0.2, 0.2, 0.2]
    @test Proper.ka_irregular_polygon_mask!(poly, xv, yv, 4, 4, 0.1f0; nsub=2) === poly
    rounded = zeros(Float32, 8, 8)
    @test Proper.ka_rounded_rectangle_mask!(rounded, 0.1f0, 0.0f0, 0.0f0, 0.1f0, 0.25f0, 0.2f0) === rounded

    table = zeros(Float32, 8, 7)
    @test Proper.ka_szoom_table!(table, 1.2f0, 8, 7, 3) === table
    tablex = zeros(Float32, 8, 7)
    tabley = zeros(Float32, 6, 7)
    @test Proper.ka_szoom_tables!(tablex, tabley, 1.2f0, 7, 3) === (tablex, tabley)
    zout = zeros(Float32, 8, 8)
    @test Proper.ka_szoom_apply!(zout, image, table, 1.2f0) === zout
    zrect = zeros(Float32, 6, 8)
    @test Proper.ka_szoom_apply!(zrect, image, tablex, tabley, 1.2f0) === zrect

    pix = zeros(Float32, 4, 4)
    @test Proper.ka_pixellate!(pix, image, 2) === pix
    @test pix[1, 1] ≈ sum(image[1:2, 1:2]) / 4
end
