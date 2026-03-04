"""Irregular polygon placeholder using bounding-box rectangle mask."""
function prop_irregular_polygon(wf::WaveFront, xverts::AbstractVector, yverts::AbstractVector; kwargs...)
    length(xverts) == length(yverts) || throw(ArgumentError("vertex arrays must have same length"))
    xmin, xmax = extrema(xverts)
    ymin, ymax = extrema(yverts)
    return prop_rectangle(wf, xmax - xmin, ymax - ymin, (xmax + xmin) / 2, (ymax + ymin) / 2; kwargs...)
end
