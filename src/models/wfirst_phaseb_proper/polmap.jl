@inline function _uniform_cubic_spline(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y size mismatch"))
    n >= 4 || throw(ArgumentError("need at least 4 samples for cubic interpolation"))

    xqf = clamp(float(xq), float(first(x)), float(last(x)))
    xf = float.(x)
    yf = float.(y)
    h = diff(xf)
    all(isapprox(hi, h[1]; rtol=0, atol=eps(h[1]) * 8) for hi in h) || throw(ArgumentError("expected uniform wavelength spacing"))
    h0 = h[1]

    A = zeros(Float64, n, n)
    b = zeros(Float64, n)

    A[1, 1] = -1 / h0
    A[1, 2] = 2 / h0
    A[1, 3] = -1 / h0
    for i in 2:(n - 1)
        him1 = h[i - 1]
        hi = h[i == n ? n - 1 : i]
        A[i, i - 1] = him1
        A[i, i] = 2 * (him1 + hi)
        A[i, i + 1] = hi
        b[i] = 6 * (((yf[i + 1] - yf[i]) / hi) - ((yf[i] - yf[i - 1]) / him1))
    end
    A[n, n - 2] = -1 / h0
    A[n, n - 1] = 2 / h0
    A[n, n] = -1 / h0

    m = A \ b
    upper = searchsortedfirst(xf, xqf)
    idx = clamp(upper - 1, 1, n - 1)
    hi = xf[idx + 1] - xf[idx]
    a = (xf[idx + 1] - xqf) / hi
    c = (xqf - xf[idx]) / hi
    return a * yf[idx] + c * yf[idx + 1] + ((a^3 - a) * m[idx] + (c^3 - c) * m[idx + 1]) * hi^2 / 6
end

function polab(polfile::AbstractString, lambda_m::Real, pupil_diam_pix::Real, condition::Integer)
    dir_out = abs(condition) == 1 ? 1 : 2
    dir_in = condition < 0 ? 1 : 2

    zamp_array = _phaseb_python_fits(polfile * "_amp.fits")
    zpha_array = _phaseb_python_fits(polfile * "_pha.fits")
    nlam = size(zamp_array, 3)
    lam_array_m = nlam == 6 ? ((0:5) .* 100 .+ 450) .* 1.0e-9 : ((0:10) .* 50 .+ 450) .* 1.0e-9
    lam = clamp(float(lambda_m), first(lam_array_m), last(lam_array_m))

    zamp = zeros(Float64, 22)
    zpha = zeros(Float64, 22)
    @inbounds for iz in 1:22
        zamp[iz] = _uniform_cubic_spline(lam_array_m, vec(Float64.(zamp_array[dir_out, dir_in, :, iz])), lam)
        zpha[iz] = _uniform_cubic_spline(lam_array_m, vec(Float64.(zpha_array[dir_out, dir_in, :, iz])), lam)
    end

    n = Int(round(float(pupil_diam_pix) * 1.1))
    n = isodd(n) ? n + 1 : n
    x = (collect(0:(n - 1)) .- n ÷ 2) ./ (float(pupil_diam_pix) / 2.0)

    amp = zeros(Float64, n, n)
    pha = zeros(Float64, n, n)

    @inbounds for j in 1:n
        y = x[j]
        r2 = x .^ 2 .+ y^2
        r = sqrt.(r2)
        r3 = r .^ 3
        r4 = r .^ 4
        r5 = r .^ 5
        r6 = r .^ 6
        t = atan.(y, x)

        for itype in 1:2
            z = itype == 1 ? zamp : zpha
            map = zeros(Float64, n)
            if itype == 1
                map .+= z[1]
            end
            map .+= z[2] .* 2 .* x
            map .+= z[3] .* 2 .* y
            map .+= z[4] .* sqrt(3) .* (2 .* r2 .- 1)
            map .+= z[5] .* sqrt(6) .* r2 .* sin.(2 .* t)
            map .+= z[6] .* sqrt(6) .* r2 .* cos.(2 .* t)
            map .+= z[7] .* sqrt(8) .* (3 .* r3 .- 2 .* r) .* sin.(t)
            map .+= z[8] .* sqrt(8) .* (3 .* r3 .- 2 .* r) .* cos.(t)
            map .+= z[9] .* sqrt(8) .* r3 .* sin.(3 .* t)
            map .+= z[10] .* sqrt(8) .* r3 .* cos.(3 .* t)
            map .+= z[11] .* sqrt(5) .* (6 .* r4 .- 6 .* r2 .+ 1)
            map .+= z[12] .* sqrt(10) .* (4 .* r4 .- 3 .* r2) .* cos.(2 .* t)
            map .+= z[13] .* sqrt(10) .* (4 .* r4 .- 3 .* r2) .* sin.(2 .* t)
            map .+= z[14] .* sqrt(10) .* r4 .* cos.(4 .* t)
            map .+= z[15] .* sqrt(10) .* r4 .* sin.(4 .* t)
            map .+= z[16] .* sqrt(12) .* (10 .* r5 .- 12 .* r3 .+ 3 .* r) .* cos.(t)
            map .+= z[17] .* sqrt(12) .* (10 .* r5 .- 12 .* r3 .+ 3 .* r) .* sin.(t)
            map .+= z[18] .* sqrt(12) .* (5 .* r5 .- 4 .* r3) .* cos.(3 .* t)
            map .+= z[19] .* sqrt(12) .* (5 .* r5 .- 4 .* r3) .* sin.(3 .* t)
            map .+= z[20] .* sqrt(12) .* r5 .* cos.(5 .* t)
            map .+= z[21] .* sqrt(12) .* r5 .* sin.(5 .* t)
            map .+= z[22] .* sqrt(7) .* (20 .* r6 .- 30 .* r4 .+ 12 .* r2 .- 1)

            if itype == 1
                amp[j, :] .= map
            else
                pha[j, :] .= map
            end
        end
    end
    return amp, pha
end

function polmap(wavefront::Proper.WaveFront, polfile::AbstractString, pupil_diam_pix::Real, condition::Integer; MUF::Real=1.0)
    n = Proper.prop_get_gridsize(wavefront)
    lambda_m = Proper.prop_get_wavelength(wavefront)

    if condition <= 2
        amp, pha = polab(polfile, lambda_m, pupil_diam_pix, condition)
    elseif condition == 5
        amp_m45_x, pha_m45_x = polab(polfile, lambda_m, pupil_diam_pix, -1)
        amp_p45_x, pha_p45_x = polab(polfile, lambda_m, pupil_diam_pix, +1)
        amp = (amp_m45_x .+ amp_p45_x) ./ 2
        pha = (pha_m45_x .+ pha_p45_x) ./ 2
    elseif condition == 6
        amp_m45_y, pha_m45_y = polab(polfile, lambda_m, pupil_diam_pix, -2)
        amp_p45_y, pha_p45_y = polab(polfile, lambda_m, pupil_diam_pix, +2)
        amp = (amp_m45_y .+ amp_p45_y) ./ 2
        pha = (pha_m45_y .+ pha_p45_y) ./ 2
    elseif condition == 10
        amp_m45_x, pha_m45_x = polab(polfile, lambda_m, pupil_diam_pix, -1)
        amp_p45_x, pha_p45_x = polab(polfile, lambda_m, pupil_diam_pix, +1)
        amp_m45_y, pha_m45_y = polab(polfile, lambda_m, pupil_diam_pix, -2)
        amp_p45_y, pha_p45_y = polab(polfile, lambda_m, pupil_diam_pix, +2)
        amp = (amp_m45_x .+ amp_p45_x .+ amp_m45_y .+ amp_p45_y) ./ 4
        pha = (pha_m45_x .+ pha_p45_x .+ pha_m45_y .+ pha_p45_y) ./ 4
    else
        throw(ArgumentError("POLMAP: unmatched condition"))
    end

    Proper.prop_multiply(wavefront, trim(amp, n))
    Proper.prop_add_phase(wavefront, trim(MUF .* pha, n))
    return wavefront
end
