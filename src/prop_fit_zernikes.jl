"""Least-squares fit of placeholder radial Zernike basis coefficients."""
function prop_fit_zernikes(wf::WaveFront, nterms::Integer)
    basis = prop_zernikes(wf, nterms)
    y = vec(real.(wf.field))
    A = reshape(basis, :, Int(nterms))
    coeff = A \ y
    return coeff
end
