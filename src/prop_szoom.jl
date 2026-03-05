@inline _szoom_round(x::T) where {T<:AbstractFloat} = x < zero(T) ? floor(x) : ceil(x)

"""Damped-sinc (Lanczos-like) zoom used by upstream PROPER `prop_magnify` default path."""
function prop_szoom(image_in::AbstractMatrix, mag0::Real, nout::Integer=0; kwargs...)
    dk = 6
    k = 13

    ny, nx = size(image_in)
    ny == nx || throw(ArgumentError("prop_szoom currently requires square input"))

    Tin = typeof(real(zero(eltype(image_in))))
    T = float(promote_type(Tin, typeof(mag0)))
    mag = T(mag0)
    mag > zero(T) || throw(ArgumentError("magnification must be positive"))

    n_in = ny
    n_out = nout == 0 ? (floor(Int, T(n_in) * mag) - k) : Int(nout)
    n_out > 0 || throw(ArgumentError("output size must be positive"))

    kk = Vector{T}(undef, k)
    @inbounds for idx in 1:k
        kk[idx] = T(idx - 1 - (k ÷ 2))
    end

    table = Matrix{T}(undef, n_out, k)
    @inbounds for i in 1:n_out
        xin = (T(i - 1 - (n_out ÷ 2))) / mag
        xphase = xin - _szoom_round(xin)
        for idx in 1:k
            x = kk[idx] - xphase
            if abs(x) <= T(dk)
                xpi = x * T(pi)
                if xpi == zero(T)
                    table[i, idx] = one(T)
                else
                    table[i, idx] = (sin(xpi) / xpi) * (sin(xpi / T(dk)) / (xpi / T(dk)))
                end
            else
                table[i, idx] = zero(T)
            end
        end
    end

    Tout = eltype(image_in) <: Complex ? Complex{T} : T
    image_out = zeros(Tout, n_out, n_out)

    @inbounds for j in 1:n_out
        yin = (T(j - 1 - (n_out ÷ 2))) / mag
        ypix = _szoom_round(yin) + T(n_in ÷ 2)
        y1 = Int(ypix) - (k ÷ 2)          # 0-based
        y2_excl = Int(ypix) + (k ÷ 2) + 1 # 0-based exclusive
        if y1 < 0 || y2_excl > n_in
            continue
        end

        for i in 1:n_out
            xin = (T(i - 1 - (n_out ÷ 2))) / mag
            xpix = _szoom_round(xin) + T(n_in ÷ 2)
            x1 = Int(xpix) - (k ÷ 2)          # 0-based
            x2_excl = Int(xpix) + (k ÷ 2) + 1 # 0-based exclusive
            if x1 < 0 || x2_excl > n_in
                continue
            end

            if eltype(image_in) <: Complex
                acc_re = zero(T)
                acc_im = zero(T)
                for co in 0:(k - 1)
                    col = x1 + co + 1
                    s_re = zero(T)
                    s_im = zero(T)
                    for ro in 0:(k - 1)
                        row = y1 + ro + 1
                        wrow = table[j, ro + 1]
                        z = image_in[row, col]
                        s_re += T(real(z)) * wrow
                        s_im += T(imag(z)) * wrow
                    end
                    wcol = table[i, co + 1]
                    acc_re += s_re * wcol
                    acc_im += s_im * wcol
                end
                image_out[j, i] = complex(acc_re, acc_im)
            else
                acc = zero(T)
                for co in 0:(k - 1)
                    col = x1 + co + 1
                    scol = zero(T)
                    for ro in 0:(k - 1)
                        row = y1 + ro + 1
                        scol += T(image_in[row, col]) * table[j, ro + 1]
                    end
                    acc += scol * table[i, co + 1]
                end
                image_out[j, i] = acc
            end
        end
    end

    return image_out
end
