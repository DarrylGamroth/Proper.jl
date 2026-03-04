"""Return placeholder mapping of Noll indices to (n,m) tuples."""
function prop_noll_zernikes(nterms::Integer)
    return [(k - 1, 0) for k in 1:Int(nterms)]
end
