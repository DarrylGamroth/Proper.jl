"""Define entrance pupil by normalizing total intensity."""
function prop_define_entrance(wf::WaveFront)
    p = sum(abs2, wf.field)
    if p > 0
        wf.field ./= sqrt(p)
    end
    return wf
end
