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

"""
    prop_begin!(field, diam, wavelength_m; beam_diam_fraction=0.5, context=nothing, workspace=nothing)

Initialize a wavefront at the entrance pupil using caller-owned complex
`field` storage.

This is the mutating counterpart to `prop_begin` for prepared or low-latency
loops that need to avoid allocating a fresh wavefront field each run.
"""
function prop_begin!(
    field::AbstractMatrix{<:Complex},
    diam::Real,
    wavelength_m::Real;
    beam_diam_fraction::Real=0.5,
    context::Union{Nothing,RunContext}=nothing,
    workspace::Union{Nothing,ProperWorkspace}=nothing,
)
    n = size(field, 1)
    n > 0 && size(field, 2) == n || throw(ArgumentError("field must be a non-empty square matrix"))
    ws = resolve_run_workspace(context, workspace)
    WT = real(eltype(field))
    λ = WT(wavelength_m)
    d = WT(diam)
    ndiam = WT(n) * WT(beam_diam_fraction)
    sampling = d / ndiam
    fill!(field, complex(one(WT), zero(WT)))
    return ws === nothing ? WaveFront(field, λ, sampling, zero(λ), d) : WaveFront(field, λ, sampling, zero(λ), d, ws)
end

"""
    prop_begin!(wf, diam, wavelength_m; beam_diam_fraction=0.5)

Reset an existing `WaveFront` in-place at the entrance pupil.

This form is intended for prepared low-latency loops that already own a
wavefront wrapper and field buffer. It avoids constructing a new `WaveFront`
object while preserving the same initialization semantics as `prop_begin!`.
"""
function prop_begin!(
    wf::WaveFront{T},
    diam::Real,
    wavelength_m::Real;
    beam_diam_fraction::Real=0.5,
) where {T<:AbstractFloat}
    field = wf.field
    n = size(field, 1)
    n > 0 && size(field, 2) == n || throw(ArgumentError("wavefront field must be a non-empty square matrix"))
    λ = T(wavelength_m)
    d = T(diam)
    ndiam = T(n) * T(beam_diam_fraction)
    sampling = d / ndiam
    fill!(field, complex(one(T), zero(T)))
    wf.wavelength_m = λ
    wf.sampling_m = sampling
    wf.z_m = zero(T)
    wf.beam_diameter_m = d
    wf.z_w0_m = zero(T)
    wf.w0_m = d / T(2)
    wf.z_rayleigh_m = pi * wf.w0_m^2 / λ
    wf.current_fratio = T(1e9)
    wf.reference_surface = PLANAR
    wf.beam_type_old = INSIDE
    wf.propagator_type = INSIDE_TO_INSIDE
    wf.rayleigh_factor = one(T)
    return wf
end
