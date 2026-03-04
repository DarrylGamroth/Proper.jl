"""Rounded-rectangle mask placeholder using plain rectangle mask."""
function prop_rounded_rectangle(wf::WaveFront, xsize::Real, ysize::Real, radius::Real=0.0, xc::Real=0.0, yc::Real=0.0; kwargs...)
    return prop_rectangle(wf, xsize, ysize, xc, yc; kwargs...)
end
