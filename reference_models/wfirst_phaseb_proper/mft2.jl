struct PhaseBMFTPlan
    left::Matrix{ComplexF64}
    right::Matrix{ComplexF64}
    scratch::Matrix{ComplexF64}
end

function prepare_phaseb_mft_plan(nfield_in::Integer, dout::Real, D::Real, nout::Integer, direction::Integer; xoffset::Real=0, yoffset::Real=0, xc::Real=0, yc::Real=0)
    x = range(-(nfield_in ÷ 2) - xc, step=1.0, length=nfield_in)
    y = range(-(nfield_in ÷ 2) - yc, step=1.0, length=nfield_in)
    du = dout / D
    u = range((-(nout ÷ 2) - xoffset / dout) * du, step=du, length=nout)
    v = range((-(nout ÷ 2) - yoffset / dout) * du, step=du, length=nout)
    xu = x * transpose(u)
    yv = y * transpose(v)
    right = (dout / D) .* exp.(-direction * 2π * im .* xu)
    left = transpose(exp.(-direction * 2π * im .* yv))
    scratch = Matrix{ComplexF64}(undef, nfield_in, nout)
    return PhaseBMFTPlan(left, right, scratch)
end

function phaseb_mft2!(out::AbstractMatrix{ComplexF64}, field_in::AbstractMatrix{<:Complex}, plan::PhaseBMFTPlan)
    mul!(plan.scratch, field_in, plan.right)
    mul!(out, plan.left, plan.scratch)
    return out
end

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
