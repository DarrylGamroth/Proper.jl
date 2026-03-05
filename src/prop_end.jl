"""Finalize propagation and return either intensity or complex field plus sampling."""
function prop_end(wf::WaveFront; noabs::Bool=false, extract::Union{Nothing,Int}=nothing)
    out = noabs ? prop_shift_center(wf.field) : prop_shift_center(abs2.(wf.field))
    if extract !== nothing
        n = extract
        ny, nx = size(out)
        n <= ny && n <= nx || throw(ArgumentError("EXTRACT exceeds array size"))
        cy = ny ÷ 2
        cx = nx ÷ 2
        r = (cy - (n ÷ 2) + 1):(cy + ((n - 1) ÷ 2) + 1)
        c = (cx - (n ÷ 2) + 1):(cx + ((n - 1) ÷ 2) + 1)
        out = @view out[r, c]
    end
    return out, wf.sampling_m
end
