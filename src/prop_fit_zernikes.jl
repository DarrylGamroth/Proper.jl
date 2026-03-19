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

    nsamples = _count_fit_zernike_samples(wavefront, pupil, pupilradius, x0, y0, RT)
    rr = Vector{RT}(undef, nsamples)
    tt = Vector{RT}(undef, nsamples)
    vv = Vector{RT}(undef, nsamples)
    _collect_fit_zernike_samples!(rr, tt, vv, wavefront, pupil, pupilradius, x0, y0)

    A = _fit_zernike_design(rr, tt, Int(nzer))
    coeff = A \ vv

    if fit
        fitmap = zeros(RT, ny, nx)
        desc = prop_noll_zernikes(Int(nzer))
        @inbounds for j in 1:nx
            x = RT(j - 1 - x0) / RT(pupilradius)
            for i in 1:ny
                y = RT(i - 1 - y0) / RT(pupilradius)
                r = hypot(x, y)
                t = atan(y, x)
                acc = zero(RT)
                for k in 1:Int(nzer)
                    acc += coeff[k] * RT(_zernike_mode(desc[k], r, t))
                end
                fitmap[i, j] = acc
            end
        end
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
