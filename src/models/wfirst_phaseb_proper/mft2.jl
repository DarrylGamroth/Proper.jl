function mft2(field_in::AbstractMatrix, dout::Real, D::Real, nout::Integer, direction::Integer; xoffset::Real=0, yoffset::Real=0, xc::Real=0, yc::Real=0)
    nfield_in = size(field_in, 2)
    x = range(-(nfield_in ÷ 2) - xc, step=1.0, length=nfield_in)
    y = range(-(nfield_in ÷ 2) - yc, step=1.0, length=nfield_in)
    du = dout / D
    u = range((-(nout ÷ 2) - xoffset / dout) * du, step=du, length=nout)
    v = range((-(nout ÷ 2) - yoffset / dout) * du, step=du, length=nout)
    xu = x * transpose(u)
    yv = y * transpose(v)
    expxu = (dout / D) .* exp.(-direction * 2π * im .* xu)
    expyv = transpose(exp.(-direction * 2π * im .* yv))
    return expyv * field_in * expxu
end
