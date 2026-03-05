"""Cubic-zoom helper compatible with PROPER's `szoom` behavior."""
function libszoom(a::AbstractMatrix, mag::Real, nout::Integer=0)
    ny, nx = size(a)
    ny == nx || throw(ArgumentError("libszoom currently requires square input"))

    out_n = nout > 0 ? Int(nout) : round(Int, ny * float(mag))
    out_n > 0 || throw(ArgumentError("nout must be positive"))
    Tin = typeof(real(zero(eltype(a))))
    T = float(promote_type(typeof(mag), Tin))

    cx_in = T(nx ÷ 2) + one(T)
    cy_in = T(ny ÷ 2) + one(T)
    cx_out = T(out_n ÷ 2)
    cy_out = T(out_n ÷ 2)

    xcoords = Vector{T}(undef, out_n)
    ycoords = Vector{T}(undef, out_n)
    invmag = inv(T(mag))
    @inbounds for j in eachindex(xcoords)
        xcoords[j] = (T(j - 1) - cx_out) * invmag + cx_in
    end
    @inbounds for i in eachindex(ycoords)
        ycoords[i] = (T(i - 1) - cy_out) * invmag + cy_in
    end

    return prop_cubic_conv(a, xcoords, ycoords; grid=true)
end
