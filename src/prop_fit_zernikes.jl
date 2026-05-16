function _fit_zernike_design(r::AbstractVector{T}, t::AbstractVector{T}, nzer::Int) where {T<:AbstractFloat}
    desc = prop_noll_zernikes(nzer)
    A = Matrix{T}(undef, length(r), nzer)
    @inbounds for j in 1:nzer
        dj = desc[j]
        for i in eachindex(r)
            A[i, j] = T(_zernike_mode(dj, r[i], t[i]))
        end
    end
    return A
end

struct FitZernikeTerm{T<:AbstractFloat}
    m::Int
    trig::UInt8
    radial_coeffs::Vector{T}
end

function _fit_zernike_terms(::Type{T}, nzer::Int) where {T<:AbstractFloat}
    descs = prop_noll_zernikes(nzer)
    terms = Vector{FitZernikeTerm{T}}(undef, nzer)
    @inbounds for k in 1:nzer
        desc = descs[k]
        top = (desc.n - desc.m) ÷ 2
        coeffs = Vector{T}(undef, top + 1)
        norm = desc.m == 0 ? sqrt(T(desc.n) + one(T)) : sqrt(T(2) * (T(desc.n) + one(T)))
        for s in 0:top
            sgn = isodd(s) ? -one(T) : one(T)
            num = sgn * T(_facti(desc.n - s))
            den = T(_facti(s)) * T(_facti((desc.n + desc.m) ÷ 2 - s)) * T(_facti((desc.n - desc.m) ÷ 2 - s))
            coeffs[s + 1] = norm * num / den
        end
        terms[k] = FitZernikeTerm{T}(desc.m, _zernike_trig_code(desc.trig), coeffs)
    end
    return terms
end

@inline function _fit_zernike_max_m(terms::AbstractVector{FitZernikeTerm{T}}) where {T<:AbstractFloat}
    max_m = 0
    @inbounds for term in terms
        max_m = max(max_m, term.m)
    end
    return max_m
end

@inline function _pow_int(x::T, n::Int) where {T<:AbstractFloat}
    acc = one(T)
    @inbounds for _ in 1:n
        acc *= x
    end
    return acc
end

@inline function _fill_zernike_angle_cache!(
    cos_m::AbstractVector{T},
    sin_m::AbstractVector{T},
    x::T,
    y::T,
    r::T,
) where {T<:AbstractFloat}
    isempty(cos_m) && return nothing
    if r == zero(T)
        fill!(cos_m, zero(T))
        fill!(sin_m, zero(T))
        return nothing
    end

    c = x / r
    s = y / r
    cm = one(T)
    sm = zero(T)
    @inbounds for m in eachindex(cos_m, sin_m)
        next_cm = cm * c - sm * s
        sm = sm * c + cm * s
        cm = next_cm
        cos_m[m] = cm
        sin_m[m] = sm
    end
    return nothing
end

@inline function _fit_zernike_radial(term::FitZernikeTerm{T}, r::T) where {T<:AbstractFloat}
    coeffs = term.radial_coeffs
    acc = coeffs[1]
    r2 = r * r
    @inbounds for k in 2:length(coeffs)
        acc = muladd(acc, r2, coeffs[k])
    end
    return term.m == 0 ? acc : acc * _pow_int(r, term.m)
end

@inline function _fit_zernike_basis_value(
    term::FitZernikeTerm{T},
    r::T,
    cos_m::AbstractVector{T},
    sin_m::AbstractVector{T},
) where {T<:AbstractFloat}
    radial = _fit_zernike_radial(term, r)
    if term.trig == _ZERNIKE_TRIG_NONE
        return radial
    elseif term.trig == _ZERNIKE_TRIG_COS
        return radial * cos_m[term.m]
    else
        return radial * sin_m[term.m]
    end
end

@inline function _fill_fit_zernike_basis!(
    basis::AbstractVector{T},
    terms::AbstractVector{FitZernikeTerm{T}},
    x::T,
    y::T,
    r::T,
    cos_m::AbstractVector{T},
    sin_m::AbstractVector{T},
) where {T<:AbstractFloat}
    _fill_zernike_angle_cache!(cos_m, sin_m, x, y, r)
    @inbounds for k in eachindex(terms, basis)
        basis[k] = _fit_zernike_basis_value(terms[k], r, cos_m, sin_m)
    end
    return basis
end

@inline function _count_fit_zernike_samples(
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius,
    x0::Integer,
    y0::Integer,
    ::Type{RT},
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    count = 0
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) / RT(pupilradius)
        for i in 1:ny
            y = RT(i - 1 - y0) / RT(pupilradius)
            r = hypot(x, y)
            count += (r < one(RT) && pupil[i, j] != 0) ? 1 : 0
        end
    end
    return count
end

@inline function _count_fit_zernike_samples(
    wavefront::AbstractMatrix,
    ::Nothing,
    pupilradius,
    x0::Integer,
    y0::Integer,
    ::Type{RT},
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    count = 0
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) / RT(pupilradius)
        for i in 1:ny
            y = RT(i - 1 - y0) / RT(pupilradius)
            count += hypot(x, y) < one(RT) ? 1 : 0
        end
    end
    return count
end

@inline function _collect_fit_zernike_samples!(
    rr::AbstractVector{RT},
    tt::AbstractVector{RT},
    vv::AbstractVector{RT},
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius,
    x0::Integer,
    y0::Integer,
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    sample = 0
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) / RT(pupilradius)
        for i in 1:ny
            y = RT(i - 1 - y0) / RT(pupilradius)
            r = hypot(x, y)
            if r < one(RT) && pupil[i, j] != 0
                sample += 1
                rr[sample] = r
                tt[sample] = atan(y, x)
                vv[sample] = RT(real(wavefront[i, j]))
            end
        end
    end
    return rr, tt, vv
end

@inline function _collect_fit_zernike_samples!(
    rr::AbstractVector{RT},
    tt::AbstractVector{RT},
    vv::AbstractVector{RT},
    wavefront::AbstractMatrix,
    ::Nothing,
    pupilradius,
    x0::Integer,
    y0::Integer,
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    sample = 0
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) / RT(pupilradius)
        for i in 1:ny
            y = RT(i - 1 - y0) / RT(pupilradius)
            r = hypot(x, y)
            if r < one(RT)
                sample += 1
                rr[sample] = r
                tt[sample] = atan(y, x)
                vv[sample] = RT(real(wavefront[i, j]))
            end
        end
    end
    return rr, tt, vv
end

function _fit_zernike_coefficients_qr(
    wavefront::AbstractMatrix,
    pupil::Union{AbstractMatrix,Nothing},
    pupilradius,
    nzer::Int,
    x0::Integer,
    y0::Integer,
    ::Type{RT},
) where {RT<:AbstractFloat}
    nsamples = _count_fit_zernike_samples(wavefront, pupil, pupilradius, x0, y0, RT)
    rr = Vector{RT}(undef, nsamples)
    tt = Vector{RT}(undef, nsamples)
    vv = Vector{RT}(undef, nsamples)
    _collect_fit_zernike_samples!(rr, tt, vv, wavefront, pupil, pupilradius, x0, y0)

    A = _fit_zernike_design(rr, tt, nzer)
    return A \ vv
end

@inline function _accumulate_fit_zernike_normal_equations!(
    gram::AbstractMatrix{RT},
    rhs::AbstractVector{RT},
    basis::AbstractVector{RT},
    terms::AbstractVector{FitZernikeTerm{RT}},
    cos_m::AbstractVector{RT},
    sin_m::AbstractVector{RT},
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius,
    x0::Integer,
    y0::Integer,
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    inv_radius = inv(RT(pupilradius))
    nzer = length(terms)
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) * inv_radius
        for i in 1:ny
            y = RT(i - 1 - y0) * inv_radius
            r = hypot(x, y)
            if r < one(RT) && pupil[i, j] != 0
                _fill_fit_zernike_basis!(basis, terms, x, y, r, cos_m, sin_m)
                v = RT(real(wavefront[i, j]))
                for a in 1:nzer
                    ba = basis[a]
                    rhs[a] += ba * v
                    for b in a:nzer
                        gram[a, b] += ba * basis[b]
                    end
                end
            end
        end
    end
    return nothing
end

@inline function _accumulate_fit_zernike_normal_equations!(
    gram::AbstractMatrix{RT},
    rhs::AbstractVector{RT},
    basis::AbstractVector{RT},
    terms::AbstractVector{FitZernikeTerm{RT}},
    cos_m::AbstractVector{RT},
    sin_m::AbstractVector{RT},
    wavefront::AbstractMatrix,
    ::Nothing,
    pupilradius,
    x0::Integer,
    y0::Integer,
) where {RT<:AbstractFloat}
    ny, nx = size(wavefront)
    inv_radius = inv(RT(pupilradius))
    nzer = length(terms)
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) * inv_radius
        for i in 1:ny
            y = RT(i - 1 - y0) * inv_radius
            r = hypot(x, y)
            if r < one(RT)
                _fill_fit_zernike_basis!(basis, terms, x, y, r, cos_m, sin_m)
                v = RT(real(wavefront[i, j]))
                for a in 1:nzer
                    ba = basis[a]
                    rhs[a] += ba * v
                    for b in a:nzer
                        gram[a, b] += ba * basis[b]
                    end
                end
            end
        end
    end
    return nothing
end

function _fit_zernike_coefficients(
    wavefront::AbstractMatrix,
    pupil::Union{AbstractMatrix,Nothing},
    pupilradius,
    nzer::Int,
    x0::Integer,
    y0::Integer,
    ::Type{RT},
) where {RT<:AbstractFloat}
    terms = _fit_zernike_terms(RT, nzer)
    max_m = _fit_zernike_max_m(terms)
    gram = zeros(RT, nzer, nzer)
    rhs = zeros(RT, nzer)
    basis = Vector{RT}(undef, nzer)
    cos_m = Vector{RT}(undef, max_m)
    sin_m = Vector{RT}(undef, max_m)

    _accumulate_fit_zernike_normal_equations!(
        gram,
        rhs,
        basis,
        terms,
        cos_m,
        sin_m,
        wavefront,
        pupil,
        pupilradius,
        x0,
        y0,
    )

    chol = cholesky!(Symmetric(gram, :U); check=false)
    if issuccess(chol)
        if cond(UpperTriangular(chol.U)) <= inv(sqrt(eps(RT)))
            coeff = ldiv!(chol, rhs)
            all(isfinite, coeff) && return coeff, terms
        end
    end

    return _fit_zernike_coefficients_qr(wavefront, pupil, pupilradius, nzer, x0, y0, RT), terms
end

function _fit_zernike_map!(
    fitmap::AbstractMatrix{RT},
    terms::AbstractVector{FitZernikeTerm{RT}},
    coeff::AbstractVector{RT},
    pupilradius,
    x0::Integer,
    y0::Integer,
) where {RT<:AbstractFloat}
    ny, nx = size(fitmap)
    inv_radius = inv(RT(pupilradius))
    max_m = _fit_zernike_max_m(terms)
    cos_m = Vector{RT}(undef, max_m)
    sin_m = Vector{RT}(undef, max_m)
    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) * inv_radius
        for i in 1:ny
            y = RT(i - 1 - y0) * inv_radius
            r = hypot(x, y)
            _fill_zernike_angle_cache!(cos_m, sin_m, x, y, r)
            acc = zero(RT)
            for k in eachindex(terms, coeff)
                acc += coeff[k] * _fit_zernike_basis_value(terms[k], r, cos_m, sin_m)
            end
            fitmap[i, j] = acc
        end
    end
    return fitmap
end

function _prop_fit_zernikes(
    ::CPUBackend,
    ::CPUBackend,
    wavefront::AbstractMatrix,
    pupil::Union{AbstractMatrix,Nothing},
    pupilradius::Real,
    nzer::Integer;
    xc::Union{Nothing,Int}=nothing,
    yc::Union{Nothing,Int}=nothing,
    fit::Bool=false,
)
    ny, nx = size(wavefront)
    pupil === nothing || size(pupil) == (ny, nx) || throw(ArgumentError("pupil size mismatch"))
    RT = pupil === nothing ?
        float(promote_type(real(eltype(wavefront)), typeof(pupilradius))) :
        float(promote_type(real(eltype(wavefront)), real(eltype(pupil)), typeof(pupilradius)))

    x0 = xc === nothing ? nx ÷ 2 : xc
    y0 = yc === nothing ? ny ÷ 2 : yc

    nterms = Int(nzer)
    coeff, terms = _fit_zernike_coefficients(wavefront, pupil, pupilradius, nterms, x0, y0, RT)

    if fit
        fitmap = Matrix{RT}(undef, ny, nx)
        _fit_zernike_map!(fitmap, terms, coeff, pupilradius, x0, y0)
        return coeff, fitmap
    end

    return coeff
end

function _prop_fit_zernikes(
    ::BackendStyle,
    ::BackendStyle,
    wavefront::AbstractMatrix,
    pupil::Union{AbstractMatrix,Nothing},
    pupilradius::Real,
    nzer::Integer;
    xc::Union{Nothing,Int}=nothing,
    yc::Union{Nothing,Int}=nothing,
    fit::Bool=false,
)
    _ = wavefront
    _ = pupil
    _ = pupilradius
    _ = nzer
    _ = xc
    _ = yc
    _ = fit
    throw(ArgumentError("prop_fit_zernikes currently supports host CPU matrices only"))
end

"""
    prop_fit_zernikes(wavefront, pupil, pupilradius, nzer; xc=nothing, yc=nothing, fit=false)

Fit Noll-ordered Zernike coefficients to a map within a pupil mask.

# Notes
- This routine currently executes on host CPU arrays.
- When `fit=true`, it also returns the reconstructed fit map.
"""
function prop_fit_zernikes(
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius::Real,
    nzer::Integer;
    xc::Union{Nothing,Int}=nothing,
    yc::Union{Nothing,Int}=nothing,
    fit::Bool=false,
)
    return _prop_fit_zernikes(
        backend_style(typeof(wavefront)),
        backend_style(typeof(pupil)),
        wavefront,
        pupil,
        pupilradius,
        nzer;
        xc=xc,
        yc=yc,
        fit=fit,
    )
end

"""Convenience fit against wavefront real part over the full-grid unit pupil."""
function prop_fit_zernikes(wf::WaveFront, nterms::Integer)
    RT = float(real(eltype(wf.field)))
    pr = RT(min(size(wf.field)...) / 2)
    return _prop_fit_zernikes(
        backend_style(typeof(wf.field)),
        CPUBackend(),
        wf.field,
        nothing,
        pr,
        Int(nterms),
    )
end
