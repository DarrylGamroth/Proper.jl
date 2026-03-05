"""Print Noll index table with (n, m, trig)."""
function prop_print_zernikes(nterms::Integer)
    for (i, d) in enumerate(prop_noll_zernikes(Int(nterms)))
        println("j=", i, " -> (n,m,trig)=(", d.n, ",", d.m, ",", d.trig, ")")
    end
    return nothing
end
