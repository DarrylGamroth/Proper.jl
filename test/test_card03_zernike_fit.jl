using Test
using Proper
using LinearAlgebra

function card03_reference_fit_zernikes_qr(
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius::Real,
    nzer::Integer;
    xc::Union{Nothing,Int}=nothing,
    yc::Union{Nothing,Int}=nothing,
    fit::Bool=false,
)
    ny, nx = size(wavefront)
    RT = float(promote_type(real(eltype(wavefront)), real(eltype(pupil)), typeof(pupilradius)))
    x0 = xc === nothing ? nx ÷ 2 : xc
    y0 = yc === nothing ? ny ÷ 2 : yc
    nterms = Int(nzer)
    coeff = Proper._fit_zernike_coefficients_qr(wavefront, pupil, pupilradius, nterms, x0, y0, RT)

    if fit
        desc = prop_noll_zernikes(nterms)
        fitmap = Matrix{RT}(undef, ny, nx)
        @inbounds for j in 1:nx
            x = RT(j - 1 - x0) / RT(pupilradius)
            for i in 1:ny
                y = RT(i - 1 - y0) / RT(pupilradius)
                r = hypot(x, y)
                t = atan(y, x)
                acc = zero(RT)
                for k in 1:nterms
                    acc += coeff[k] * RT(Proper._zernike_mode(desc[k], r, t))
                end
                fitmap[i, j] = acc
            end
        end
        return coeff, fitmap
    end

    return coeff
end

function card03_circular_mask(n::Int, radius::Float64, x0::Int, y0::Int; clipped::Bool=false)
    mask = zeros(Float64, n, n)
    @inbounds for j in 1:n
        x = (j - 1 - x0) / radius
        for i in 1:n
            y = (i - 1 - y0) / radius
            missing = clipped && i > y0 && j < x0 - n ÷ 10
            mask[i, j] = hypot(x, y) < 1.0 && !missing ? 1.0 : 0.0
        end
    end
    return mask
end

function card03_wavefront_map(n::Int, nterms::Int)
    wf = prop_begin(1.0, 550e-9, n)
    coeffs = [((isodd(k) ? 1.0 : -1.0) * (0.5 + k / nterms) * 1e-9) for k in 1:nterms]
    zmap = prop_zernikes(wf, collect(1:nterms), coeffs; no_apply=true)
    @inbounds for j in 1:n
        for i in 1:n
            zmap[i, j] += 2e-10 * sinpi(i / n) * cospi(2j / n)
        end
    end
    return zmap
end

@testset "Card 03 Zernike fit active-sample path" begin
    cases = (
        (grid=32, nzer=8, radius=14.0, xc=16, yc=16, clipped=false),
        (grid=64, nzer=22, radius=28.0, xc=33, yc=30, clipped=true),
    )

    for case in cases
        wavefront = card03_wavefront_map(case.grid, max(case.nzer, 22))
        pupil = card03_circular_mask(case.grid, case.radius, case.xc, case.yc; clipped=case.clipped)

        coeff_ref = card03_reference_fit_zernikes_qr(
            wavefront,
            pupil,
            case.radius,
            case.nzer;
            xc=case.xc,
            yc=case.yc,
        )
        coeff = prop_fit_zernikes(
            wavefront,
            pupil,
            case.radius,
            case.nzer;
            xc=case.xc,
            yc=case.yc,
        )
        @test coeff ≈ coeff_ref rtol = 1e-9 atol = 5e-20

        coeff_fit_ref, fitmap_ref = card03_reference_fit_zernikes_qr(
            wavefront,
            pupil,
            case.radius,
            case.nzer;
            xc=case.xc,
            yc=case.yc,
            fit=true,
        )
        coeff_fit, fitmap = prop_fit_zernikes(
            wavefront,
            pupil,
            case.radius,
            case.nzer;
            xc=case.xc,
            yc=case.yc,
            fit=true,
        )
        @test coeff_fit ≈ coeff_fit_ref rtol = 1e-9 atol = 5e-20
        @test fitmap ≈ fitmap_ref rtol = 1e-9 atol = 5e-20
        @test size(fitmap) == size(wavefront)
        @test fitmap isa Matrix{Float64}
    end

    sparse_wavefront = card03_wavefront_map(16, 8)
    sparse_pupil = zeros(Float64, 16, 16)
    sparse_pupil[8, 8] = 1.0
    sparse_ref = card03_reference_fit_zernikes_qr(sparse_wavefront, sparse_pupil, 7.0, 4; xc=8, yc=8)
    sparse_coeff = prop_fit_zernikes(sparse_wavefront, sparse_pupil, 7.0, 4; xc=8, yc=8)
    @test sparse_coeff ≈ sparse_ref rtol = 1e-9 atol = 5e-20
end
