"""Initialize a wavefront at entrance pupil."""
function prop_begin(
    diam::Real,
    wavelength_m::Real,
    gridsize::Integer;
    beam_diam_fraction::Real=0.5,
    context::Union{Nothing,RunContext}=nothing,
    workspace::Union{Nothing,ProperWorkspace}=nothing,
)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    ws = resolve_run_workspace(context, workspace)
    WT = something(workspace_float_type(ws), float(typeof(wavelength_m)))
    λ = WT(wavelength_m)
    d = WT(diam)
    ndiam = WT(n) * WT(beam_diam_fraction)
    sampling = d / ndiam
    field = unit_field(WT, n, ws)
    return ws === nothing ? WaveFront(field, λ, sampling, zero(λ), d) : WaveFront(field, λ, sampling, zero(λ), d, ws)
end
