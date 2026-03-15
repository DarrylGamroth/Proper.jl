function mft2(field_in::AbstractMatrix, dout::Real, D::Real, nout::Integer, direction::Integer; xoffset::Real=0, yoffset::Real=0, xc::Real=0, yc::Real=0)
    nfield_in = size(field_in, 2)
    x = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2) .- xc)
    y = (collect(0:(nfield_in - 1)) .- (nfield_in ÷ 2) .- yc)
    u = (collect(0:(nout - 1)) .- (nout ÷ 2) .- xoffset / dout) .* (dout / D)
    v = (collect(0:(nout - 1)) .- (nout ÷ 2) .- yoffset / dout) .* (dout / D)
    xu = x * transpose(u)
    yv = y * transpose(v)
    expxu = (dout / D) .* exp.(-direction * 2π * im .* xu)
    expyv = transpose(exp.(-direction * 2π * im .* yv))
    return expyv * field_in * expxu
end
