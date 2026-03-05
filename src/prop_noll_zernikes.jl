"""Return Noll-ordered Zernike descriptors as `(n, m, trig)` tuples."""
function prop_noll_zernikes(nterms::Integer)
    nterms > 0 || throw(ArgumentError("nterms must be positive"))
    out = Vector{NamedTuple{(:n, :m, :trig),Tuple{Int,Int,Symbol}}}(undef, nterms)

    iz = 1
    n = 0
    while iz <= nterms
        for m in (n % 2):2:n
            if n == 0
                out[iz] = (n=0, m=0, trig=:none)
                iz += 1
                iz > nterms && break
                continue
            end

            if m == 0
                out[iz] = (n=n, m=0, trig=:none)
                iz += 1
                iz > nterms && break
            else
                out[iz] = (n=n, m=m, trig=iseven(iz) ? :cos : :sin)
                iz += 1
                iz > nterms && break
                out[iz] = (n=n, m=m, trig=iseven(iz) ? :cos : :sin)
                iz += 1
                iz > nterms && break
            end
        end
        n += 1
    end

    return out
end
