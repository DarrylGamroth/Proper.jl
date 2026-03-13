const SZOOM_DK = 6
const SZOOM_K = 13

@inline _szoom_round(x::T) where {T<:AbstractFloat} = x < zero(T) ? floor(x) : ceil(x)

@inline function _szoom_magtype(image_in::AbstractMatrix, mag0::Real)
    Tin = typeof(real(zero(eltype(image_in))))
    return float(promote_type(Tin, typeof(mag0)))
end

@inline function _szoom_output_eltype(image_in::AbstractMatrix, ::Type{T}) where {T<:AbstractFloat}
    return eltype(image_in) <: Complex ? Complex{T} : T
end

@inline function _szoom_output_size(n_in::Integer, mag::T, nout::Integer) where {T<:AbstractFloat}
    return nout == 0 ? (floor(Int, T(n_in) * mag) - SZOOM_K) : Int(nout)
end

@inline function _fill_szoom_table_loop!(table::AbstractMatrix{T}, mag::T) where {T<:AbstractFloat}
    n_out, k = size(table)
    @inbounds for i in 1:n_out
        xin = T(i - 1 - (n_out ÷ 2)) / mag
        xphase = xin - _szoom_round(xin)
        for idx in 1:k
            x = T(idx - 1 - (k ÷ 2)) - xphase
            if abs(x) <= T(SZOOM_DK)
                xpi = x * T(pi)
                if xpi == zero(T)
                    table[i, idx] = one(T)
                else
                    table[i, idx] = (sin(xpi) / xpi) * (sin(xpi / T(SZOOM_DK)) / (xpi / T(SZOOM_DK)))
                end
            else
                table[i, idx] = zero(T)
            end
        end
    end
    return table
end

function _apply_szoom_loop!(out::AbstractMatrix, image_in::AbstractMatrix, table::AbstractMatrix{T}, mag::T) where {T<:AbstractFloat}
    n_in = size(image_in, 1)
    n_out = size(out, 1)
    fill!(out, zero(eltype(out)))

    @inbounds for j in 1:n_out
        yin = T(j - 1 - (n_out ÷ 2)) / mag
        ypix = _szoom_round(yin) + T(n_in ÷ 2)
        y1 = Int(ypix) - (SZOOM_K ÷ 2)
        y2_excl = Int(ypix) + (SZOOM_K ÷ 2) + 1
        if y1 < 0 || y2_excl > n_in
            continue
        end

        for i in 1:n_out
            xin = T(i - 1 - (n_out ÷ 2)) / mag
            xpix = _szoom_round(xin) + T(n_in ÷ 2)
            x1 = Int(xpix) - (SZOOM_K ÷ 2)
            x2_excl = Int(xpix) + (SZOOM_K ÷ 2) + 1
            if x1 < 0 || x2_excl > n_in
                continue
            end

            if eltype(image_in) <: Complex
                acc_re = zero(T)
                acc_im = zero(T)
                for co in 0:(SZOOM_K - 1)
                    col = x1 + co + 1
                    s_re = zero(T)
                    s_im = zero(T)
                    for ro in 0:(SZOOM_K - 1)
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
                out[j, i] = complex(acc_re, acc_im)
            else
                acc = zero(T)
                for co in 0:(SZOOM_K - 1)
                    col = x1 + co + 1
                    scol = zero(T)
                    for ro in 0:(SZOOM_K - 1)
                        row = y1 + ro + 1
                        scol += T(image_in[row, col]) * table[j, ro + 1]
                    end
                    acc += scol * table[i, co + 1]
                end
                out[j, i] = acc
            end
        end
    end

    return out
end

function _prop_szoom!(
    ::SamplingLoopExecStyle,
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag::T,
) where {T<:AbstractFloat}
    table = Matrix{T}(undef, size(out, 1), SZOOM_K)
    _fill_szoom_table_loop!(table, mag)
    return _apply_szoom_loop!(out, image_in, table, mag)
end

function _prop_szoom!(
    ::SamplingKAExecStyle,
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag::T,
) where {T<:AbstractFloat}
    table = similar(out, T, size(out, 1), SZOOM_K)
    ka_szoom_table!(table, mag, size(out, 1), SZOOM_K, SZOOM_DK)
    return ka_szoom_apply!(out, image_in, table, mag)
end

function _prop_szoom!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag::T,
) where {T<:AbstractFloat}
    ny, nx = size(image_in)
    ny == nx || throw(ArgumentError("prop_szoom currently requires square input"))
    size(out, 1) == size(out, 2) || throw(ArgumentError("output must be square"))
    return _prop_szoom!(sampling_exec_style(typeof(out), size(out, 1), size(out, 2)), out, image_in, mag)
end

function prop_szoom!(
    out::AbstractMatrix,
    image_in::AbstractMatrix,
    mag0::Real,
)
    T = _szoom_magtype(image_in, mag0)
    mag = T(mag0)
    mag > zero(T) || throw(ArgumentError("magnification must be positive"))
    return _prop_szoom!(out, image_in, mag)
end

"""Damped-sinc (Lanczos-like) zoom used by upstream PROPER `prop_magnify` default path."""
function prop_szoom(image_in::AbstractMatrix, mag0::Real, nout::Integer=0; kwargs...)
    T = _szoom_magtype(image_in, mag0)
    mag = T(mag0)
    mag > zero(T) || throw(ArgumentError("magnification must be positive"))

    ny, nx = size(image_in)
    ny == nx || throw(ArgumentError("prop_szoom currently requires square input"))

    n_out = _szoom_output_size(ny, mag, nout)
    n_out > 0 || throw(ArgumentError("output size must be positive"))

    Tout = _szoom_output_eltype(image_in, T)
    out = similar(image_in, Tout, n_out, n_out)
    return _prop_szoom!(out, image_in, mag)
end
