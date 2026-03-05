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

"""Fit Noll-ordered Zernike coefficients to a map within pupil mask."""
function prop_fit_zernikes(
    wavefront::AbstractMatrix,
    pupil::AbstractMatrix,
    pupilradius::Real,
    nzer::Integer;
    xc::Union{Nothing,Int}=nothing,
    yc::Union{Nothing,Int}=nothing,
    fit::Bool=false,
)
    ny, nx = size(wavefront)
    size(pupil) == (ny, nx) || throw(ArgumentError("pupil size mismatch"))
    RT = float(promote_type(real(eltype(wavefront)), real(eltype(pupil)), typeof(pupilradius)))

    x0 = xc === nothing ? nx ÷ 2 : xc
    y0 = yc === nothing ? ny ÷ 2 : yc

    rr = Vector{RT}()
    tt = Vector{RT}()
    vv = Vector{RT}()

    @inbounds for j in 1:nx
        x = RT(j - 1 - x0) / RT(pupilradius)
        for i in 1:ny
            y = RT(i - 1 - y0) / RT(pupilradius)
            r = hypot(x, y)
            if r < one(RT) && pupil[i, j] != 0
                push!(rr, r)
                push!(tt, atan(y, x))
                push!(vv, RT(real(wavefront[i, j])))
            end
        end
    end

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

"""Convenience fit against wavefront real part over full grid pupil."""
function prop_fit_zernikes(wf::WaveFront, nterms::Integer)
    ny, nx = size(wf.field)
    RT = float(real(eltype(wf.field)))
    pupil = ones(RT, ny, nx)
    pr = RT(min(ny, nx) / 2)
    return prop_fit_zernikes(real.(wf.field), pupil, pr, Int(nterms))
end
