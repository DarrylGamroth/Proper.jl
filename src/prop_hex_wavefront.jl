"""Apply segmented hex aperture transmission and optional per-segment hex-Zernike phase."""
function prop_hex_wavefront(
    wf::WaveFront,
    nrings::Integer,
    hexrad::Real,
    hexsep::Real,
    zernike_val=0;
    xcenter::Real=0.0,
    ycenter::Real=0.0,
    darkcenter::Bool=false,
    no_apply::Bool=false,
    omit=Int[],
    rotation::Real=0.0,
    kwargs...,
)
    xc = haskey(kwargs, :XCENTER) ? float(kwargs[:XCENTER]) : float(xcenter)
    yc = haskey(kwargs, :YCENTER) ? float(kwargs[:YCENTER]) : float(ycenter)
    dark = darkcenter || switch_set(:DARKCENTER; kwargs...)
    noapply = no_apply || switch_set(:NO_APPLY; kwargs...)
    angle = haskey(kwargs, :ROTATION) ? float(kwargs[:ROTATION]) : float(rotation)
    omit_ids = haskey(kwargs, :OMIT) ? Int.(collect(kwargs[:OMIT])) : Int.(collect(omit))
    omit_set = Set(omit_ids)

    n = prop_get_gridsize(wf)
    RT = real(eltype(wf.field))
    aperture = zeros(RT, n, n)

    include_phase = !(zernike_val isa Number && iszero(zernike_val))
    phase = zeros(RT, n, n)
    totseg = zeros(RT, n, n)

    zvals = nothing
    if include_phase
        zraw = float.(collect(zernike_val))
        nhex = nrings * (nrings + 1) * 3 + 1
        if ndims(zraw) != 2
            throw(ArgumentError("zernike_val must be a 2-D array with 22 coefficients per segment"))
        end
        if size(zraw, 1) == 22 && size(zraw, 2) == nhex
            zvals = zraw
        elseif size(zraw, 1) == nhex && size(zraw, 2) == 22
            zvals = permutedims(zraw)
        else
            throw(ArgumentError("zernike_val must have shape (22,nhex) or (nhex,22), nhex=$(nhex)"))
        end
    end

    segi = 0
    ang = deg2rad(angle)
    s, c = sincos(ang)

    @inbounds for iring in 0:nrings
        x = hexsep * cospi(1 / 6) * iring
        y = -nrings * hexsep + iring * hexsep * 0.5

        for iseg in 0:(2 * nrings - iring)
            if !(segi in omit_set)
                xhex = x * c - y * s + xc
                yhex = x * s + y * c + yc

                segment = prop_polygon(wf, 6, hexrad, xhex, yhex; ROTATION=angle)
                if !(iring == 0 && iseg == nrings && dark)
                    aperture .+= segment
                end

                if include_phase
                    segmask = segment .!= 0
                    coeff = @view(zvals[:, segi + 1])
                    phase .+= segmask .* prop_hex_zernikes(1:22, coeff, n, prop_get_sampling(wf), hexrad, xhex, yhex; rotation=angle)
                    totseg .+= segmask
                end
            end
            segi += 1

            if iring != 0
                if !(segi in omit_set)
                    xhex = -x * c - y * s + xc
                    yhex = -x * s + y * c + yc

                    segment = prop_polygon(wf, 6, hexrad, xhex, yhex; ROTATION=angle)
                    aperture .+= segment

                    if include_phase
                        segmask = segment .!= 0
                        coeff = @view(zvals[:, segi + 1])
                        phase .+= segmask .* prop_hex_zernikes(1:22, coeff, n, prop_get_sampling(wf), hexrad, xhex, yhex; rotation=angle)
                        totseg .+= segmask
                    end
                end
                segi += 1
            end

            y += hexsep
        end
    end

    if include_phase
        @inbounds for j in axes(phase, 2)
            for i in axes(phase, 1)
                t = totseg[i, j]
                if t != 0
                    phase[i, j] /= t
                end
            end
        end
    end

    if !noapply
        prop_multiply(wf, aperture)
        if include_phase
            prop_add_phase(wf, phase)
        end
    end

    return include_phase ? (aperture, phase) : aperture
end
