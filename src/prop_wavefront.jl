"""Create a wavefront container with a complex field initialized to 1+0im."""
function prop_wavefront(
    gridsize::Integer,
    wavelength_m::Real,
    beam_diameter_m::Real;
    sampling_m::Union{Nothing,Real}=nothing,
    context::Union{Nothing,RunContext}=nothing,
    workspace::Union{Nothing,ProperWorkspace}=nothing,
)
    n = Int(gridsize)
    n > 0 || throw(ArgumentError("gridsize must be positive"))
    ws = resolve_run_workspace(context, workspace)
    WT = something(workspace_float_type(ws), float(typeof(wavelength_m)))
    λ = WT(wavelength_m)
    d = WT(beam_diameter_m)
    s = sampling_m === nothing ? d / WT(n) : WT(sampling_m)
    field = unit_field(WT, n, ws)
    return ws === nothing ? WaveFront(field, λ, s, zero(λ), d) : WaveFront(field, λ, s, zero(λ), d, ws)
end
