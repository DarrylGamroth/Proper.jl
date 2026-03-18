using Proper

const CORE_PROPAGATION_TAIL_GRID_N = 512
const CORE_PROPAGATION_TAIL_SAMPLES = 20
const CORE_PROPAGATION_TAIL_PROFILE_ITERS = 8

const CORE_PROPAGATION_TAIL_SEQUENCE = (
    ("lens_1", :lens, 1.25),
    ("prop_1", :propagate, 0.18),
    ("prop_2", :propagate, 0.22),
    ("lens_2", :lens, 0.90),
    ("prop_3", :propagate, 0.14),
    ("prop_4", :propagate, 0.31),
    ("lens_3", :lens, 0.75),
    ("prop_5", :propagate, 0.27),
)

"""
    core_propagation_tail!(wf, ctx)

Run a small synthetic propagation sequence dominated by repeated
`prop_lens`/`prop_propagate` transitions. This workload exists to benchmark and
profile shared core propagation behavior without using the WFIRST model itself
as the optimization target.
"""
function core_propagation_tail!(wf::WaveFront, ctx::RunContext)
    @inbounds for (_, op, value) in CORE_PROPAGATION_TAIL_SEQUENCE
        if op === :lens
            prop_lens(wf, value, ctx)
        else
            prop_propagate(wf, value, ctx)
        end
    end
    return wf
end

function cpu_core_propagation_tail_case(::Type{T}=Float64, grid_n::Integer=CORE_PROPAGATION_TAIL_GRID_N) where {T<:AbstractFloat}
    wf = prop_begin(T(2.4), T(0.55e-6), grid_n; beam_diam_fraction=T(0.5))
    prop_circular_aperture(wf, T(0.6))
    prop_define_entrance(wf)
    ctx = RunContext(wf)
    snap = capture_wavefront_state(wf)
    return wf, ctx, snap
end

function gpu_core_propagation_tail_case(
    builder::F,
    sync!::S,
    ::Type{T}=Float64,
    grid_n::Integer=CORE_PROPAGATION_TAIL_GRID_N,
) where {F<:Function,S<:Function,T<:AbstractFloat}
    wf = builder(T, 2.4, 0.55e-6, grid_n; beam_diam_fraction=T(0.5))
    ctx = RunContext(wf)
    prop_circular_aperture(wf, T(0.6))
    prop_define_entrance(wf)
    sync!()
    snap = capture_wavefront_state(wf)
    return wf, ctx, snap
end

function restore_and_run_core_propagation_tail!(wf::WaveFront, snap, ctx::RunContext)
    restore_wavefront_state!(wf, snap)
    return core_propagation_tail!(wf, ctx)
end
