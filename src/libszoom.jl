"""Cubic-zoom helper compatible with PROPER's `szoom` behavior."""
function libszoom(a::AbstractMatrix, mag::Real, nout::Integer=0)
    ny, nx = size(a)
    ny == nx || throw(ArgumentError("libszoom currently requires square input"))

    out_n = nout > 0 ? Int(nout) : round(Int, ny * float(mag))
    cy_in = (ny + 1) / 2
    cx_in = (nx + 1) / 2
    cy_out = (out_n + 1) / 2
    cx_out = (out_n + 1) / 2

    xcoords = Vector{Float64}(undef, out_n)
    ycoords = Vector{Float64}(undef, out_n)
    @inbounds for j in eachindex(xcoords)
        xcoords[j] = (j - cx_out) / mag + cx_in
    end
    @inbounds for i in eachindex(ycoords)
        ycoords[i] = (i - cy_out) / mag + cy_in
    end

    return prop_cubic_conv(a, xcoords, ycoords; grid=true)
end
