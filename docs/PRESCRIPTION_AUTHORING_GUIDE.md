# Prescription Authoring Guide

## Purpose
This guide is for users who want to:

- port an existing Python or MATLAB PROPER prescription to Julia
- write a new prescription directly against `Proper.jl`
- decide when to stay on the familiar `prop_*` surface and when to use prepared execution

It complements:

- `docs/MIGRATION_GUIDE.md` for cross-language compatibility expectations
- `docs/PREPARED_EXECUTION_GUIDE.md` for repeated-run execution patterns
- `docs/api_contract.md` for the stable public API contract

## Recommended Starting Point
Start with the familiar PROPER shape:

```julia
using Proper

function my_prescription(λm, n; kwargs...)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)
    prop_lens(wf, 10.0)
    prop_propagate(wf, 10.0)
    return prop_end(wf)
end
```

Then run it with:

```julia
psf, sampling = prop_run(my_prescription, 0.55, 256)
```

That is the right starting point even if you expect to move to prepared execution
later. It keeps parity validation simple and makes upstream comparisons easier.

## Authoring Model

### Public Boundary
At the public boundary, keep the familiar PROPER conventions:

- wavelength passed to `prop_run` in microns
- grid size passed as an integer dimension
- prescription returns either:
  - `prop_end(wf)` output, or
  - an explicit `(psf, sampling)` tuple

### Inside The Prescription
Inside a prescription:

- work with the `WaveFront`
- use the public `prop_*` routines
- keep the flow explicit and readable
- prefer ordinary Julia control flow over dynamic string/symbol dispatch

That means a Julia prescription should usually look like a normal Julia function,
not like a transliterated Python module body.

## Porting From Python Or MATLAB

### What Usually Ports Directly
Most optical sequence code ports with minimal structural change:

- `prop_begin`
- apertures/obscurations
- `prop_define_entrance`
- `prop_lens`
- `prop_propagate`
- `prop_errormap`
- `prop_dm`
- `prop_end`

The biggest changes are usually around:

- Python module globals
- `PASSVALUE` shape and keyword handling
- plotting
- FITS I/O
- repeated-run execution structure

### What Not To Copy Literally
Do not preserve upstream implementation quirks just because they exist in the
source tree.

In particular:

- do not recreate Python module-global mutable state
- do not keep runtime `eval`/dynamic import patterns inside prescriptions
- do not branch on string/symbol mode flags in hot loops
- do not transliterate line-by-line when Julia has a clearer direct expression

Parity-first does not mean transliteration-first.

### Keyword Translation
Compatibility keywords remain accepted in uppercase form, but Julia code should
prefer ordinary keyword arguments where practical.

Typical translation pattern:

```julia
function my_prescription(λm, n; PASSVALUE=nothing, use_dm=false)
    wf = prop_begin(1.0, λm, n)
    if use_dm
        # ...
    end
    return prop_end(wf)
end
```

Use `PASSVALUE` when you need upstream-style compatibility or parity harness
alignment. Prefer explicit Julia keywords for new code.

## `PASSVALUE` Guidance
`PASSVALUE` is useful as a compatibility adapter, not as the ideal long-term
shape for all Julia code.

Use it when:

- porting an upstream prescription
- matching an existing parity harness
- carrying a heterogeneous set of optional inputs through `prop_run`

Prefer explicit keywords or prepared assets when:

- the model is now Julia-native
- the same assets are reused across runs
- the call surface is stable enough to type explicitly

## Use Prepared Execution Deliberately
Do not start with prepared execution just because it exists.

Use the execution surfaces in this order:

1. `prop_run(...)`
2. `prepare_prescription(...)`
3. `prepare_prescription_batch(...)`
4. `prepare_model(...)` and `prepare_asset_pool(...)`

That ordering matters because:

- parity validation is simplest with plain `prop_run`
- prepared execution is an ownership/performance tool
- prepared execution should not hide prescription logic during the initial port

## Prescription Structure Recommendations

### Keep Optical Steps Obvious
Prefer this:

```julia
wf = prop_begin(beam_diameter, λm, n)
prop_circular_aperture(wf, radius)
prop_define_entrance(wf)
prop_lens(wf, fl_oap1)
prop_propagate(wf, d_oap1_oap2)
return prop_end(wf)
```

Over collapsing the whole prescription into helpers too early. The optical
sequence should remain easy to audit against upstream references.

### Push Reuse Into Helpers, Not Hidden State
If several prescriptions share setup logic, extract normal Julia helper
functions. Keep them explicit:

```julia
function apply_common_front_end!(wf)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)
    return wf
end
```

That is preferable to hidden mutable module state.

### Use Mutating APIs Intentionally
For prescription authoring, the non-mutating public wrappers are usually the
right starting point.

Reach for `...!` variants when:

- you own the destination array
- profiling shows the allocation matters
- the mutating form keeps the prescription readable

Do not force mutating style into every prescription just for appearances.

## Arrays, Backends, And GPU Readiness
When writing reusable prescription code:

- accept normal Julia values and arrays at the boundary
- avoid forcing `Array` unless host materialization is actually required
- avoid scalar indexing assumptions that break GPU execution
- keep array-manipulation helpers generic over `AbstractArray` when practical

If a prescription needs backend-specific assets, isolate that in preparation
logic rather than scattering backend checks through the optical sequence.

## FITS And External Data
For map-heavy prescriptions:

- load assets through the established FITS helpers
- keep array-order assumptions explicit
- validate against the Python baseline when parity matters
- use MATLAB as the semantic reference for column-major correctness questions

If a model requires compatibility reshaping or metadata normalization, keep that
policy local to the model layer and document it.

## Testing A Prescription
For a new or ported prescription, validate in this order:

1. smoke test the Julia prescription directly
2. compare against a Python parity case if an upstream prescription exists
3. move repeated runs onto prepared execution if needed
4. benchmark only after parity is established

Useful checks:

- return shape is stable
- sampling matches expected units
- complex-field vs intensity output is intentional
- random paths are seeded in tests

## Common Mistakes
- porting upstream globals directly instead of using function arguments
- treating `PASSVALUE` as the only valid Julia interface
- materializing host arrays inside otherwise generic code
- adding model-specific optimization tricks before parity is established
- rewriting a prescription around prepared execution before the plain `prop_run`
  path is validated

## Minimal Migration Pattern
For an existing upstream prescription, the usual workflow is:

1. copy the file structure and prescription name for traceability
2. port the optical sequence using familiar `prop_*` calls
3. get `prop_run(...)` working first
4. confirm parity against the Python baseline
5. only then introduce prepared execution or cached assets if repeated runs need it

## Related Docs
- `README.md`
- `docs/MIGRATION_GUIDE.md`
- `docs/PREPARED_EXECUTION_GUIDE.md`
- `docs/API_EXAMPLES.md`
- `docs/api_contract.md`
- `docs/compat_decisions.md`
