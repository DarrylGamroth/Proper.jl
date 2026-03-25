# Migration Guide: Python/MATLAB PROPER to Proper.jl

## Read This First
If you already know upstream PROPER, this is the one document you should read
first.

The intended migration path is:

1. port the prescription with the familiar `prop_*` calls
2. run it with `prop_run(...)`
3. keep `PASSVALUE` only when you need upstream-style compatibility
4. introduce prepared execution only after the plain port is correct

Supporting documents after this one:
- [Runnable API examples](API_EXAMPLES.md)
- [Prepared execution guide](PREPARED_EXECUTION_GUIDE.md)
- [Prescription authoring guide](PRESCRIPTION_AUTHORING_GUIDE.md)

## Stable Expectations
- Public routine naming remains familiar (`prop_*`).
- Python 3.3.4 remains the executable parity baseline.
- MATLAB 3.3.1 remains the semantic reference.
- FITS IO uses `FITSIO.jl`.
- The public execution contract is still tuple-based:
  - `prop_run(...) -> (psf, sampling)`
  - `prop_run_multi(...) -> (stack, samplings)`

## Upstream PROPER Cheat Sheet

| Upstream concept | Proper.jl equivalent | When to use it |
| --- | --- | --- |
| `prop_run(prescription, λ, n, PASSVALUE=...)` | `prop_run(prescription, λ, n; PASSVALUE=...)` | default migration path |
| repeated single prescription runs | `prepare_prescription(...)` | reuse one normalized call shape |
| repeated or parallel runs with varying `PASSVALUE` | `prepare_prescription_batch(...)` | slot-local context reuse |
| one named reusable configured model | `prepare_model(...)` | application-facing execution object |
| wavelength sweep | `prop_run_multi(runs)` where `runs` is a vector of prepared runs | fixed-shape throughput path |
| uppercase compatibility keywords | still accepted | migration and parity |
| lowercase Julia keywords | also accepted | Julia-native code |

The practical rule is:
- start exactly where upstream usage starts
- only introduce prepared execution when repeated work or cached state is the
  actual problem

## Prescription Signature Translation

The most common upstream prescription shapes map to Julia like this.

### Plain prescription

```julia
function my_prescription(λm, n)
    wf = prop_begin(1.0, λm, n)
    return prop_end(wf)
end
```

### Upstream-style `PASSVALUE`

```julia
function my_prescription(λm, n; PASSVALUE=nothing)
    wf = prop_begin(1.0, λm, n)
    return prop_end(wf)
end
```

### Positional pass value

`prop_run` also accepts the common upstream calling shape where the
prescription takes a positional third argument:

```julia
function my_prescription(λm, n, passvalue)
    wf = prop_begin(1.0, λm, n)
    return prop_end(wf)
end
```

`prop_run` normalizes both forms at the API boundary.

## Side-By-Side Example

### Python

```python
def simple_prescription(lam, n, PASSVALUE=None):
    wf = proper.prop_begin(1.0, lam, n)
    proper.prop_circular_aperture(wf, 0.5)
    proper.prop_define_entrance(wf)
    proper.prop_lens(wf, 10.0)
    proper.prop_propagate(wf, 10.0)
    return proper.prop_end(wf)
```

### MATLAB

```matlab
function [psf, sampling] = simple_prescription(lam, n, optval)
    wf = prop_begin(1.0, lam, n);
    wf = prop_circular_aperture(wf, 0.5);
    wf = prop_define_entrance(wf);
    wf = prop_lens(wf, 10.0);
    wf = prop_propagate(wf, 10.0);
    [psf, sampling] = prop_end(wf);
end
```

### Julia

```julia
using Proper

function simple_prescription(λm, n; PASSVALUE=nothing)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)
    prop_lens(wf, 10.0)
    prop_propagate(wf, 10.0)
    return prop_end(wf)
end

psf, sampling = prop_run(simple_prescription, 0.55, 256)
```

What changes:
- Julia uses normal function definitions instead of Python module-level state.
- `PASSVALUE` is a keyword by convention, though positional third-argument
  prescriptions also work.
- The optical sequence stays almost identical.

What does not change:
- wavelength at the public boundary is still in microns
- the prescription still returns `prop_end(wf)` or `(psf, sampling)`
- the familiar `prop_*` surface remains the main way to write the prescription

## DM / FITS / Map-Driven Prescription Pattern

This is the common nontrivial migration shape many upstream users care about.

```julia
using Proper

function dm_map_prescription(λm, n; PASSVALUE=nothing)
    pass = PASSVALUE === nothing ? NamedTuple() : PASSVALUE

    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)

    if haskey(pass, :errormap_path)
        prop_errormap(wf, pass[:errormap_path], WAVEFRONT=true)
    end

    if haskey(pass, :dm_map)
        prop_dm(wf, pass[:dm_map])
    end

    prop_lens(wf, 10.0)
    prop_propagate(wf, 10.0)
    return prop_end(wf)
end
```

Run it with:

```julia
psf, sampling = prop_run(
    dm_map_prescription,
    0.55,
    256;
    PASSVALUE=Dict(
        :errormap_path => "phase_map.fits",
        :dm_map => zeros(256, 256),
    ),
)
```

Important behavior:
- FITS files are read on the host, then promoted to the wavefront backend where
  feasible.
- `prop_dm(wf, dm_map)` expects a map already sampled on the wavefront grid.
- Use `prop_dm` with the full upstream-compatible actuator-space form only when
  you are intentionally modeling actuator geometry and influence functions.

## Keyword Translation Rules

- Uppercase compatibility keywords are accepted.
- Lowercase Julia aliases are also accepted at the public boundary.
- If you are writing new Julia-native code, prefer ordinary Julia keywords.

Examples:

```julia
rot = prop_rotate(img, 15.0; METH="cubic")
rot = prop_rotate(img, 15.0; meth="cubic")
```

```julia
mag = prop_magnify(img, 1.5; QUICK=true)
mag = prop_magnify(img, 1.5; quick=true)
```

## `PASSVALUE` Guidance

Use `PASSVALUE` when:
- porting an upstream prescription
- matching an existing parity harness
- carrying a heterogeneous set of optional inputs through `prop_run`

Prefer explicit Julia keywords or prepared assets when:
- the model is now Julia-native
- the same assets are reused across runs
- the call surface is stable enough to type explicitly

## Prepared Execution Mental Model

Prepared execution is not a compatibility mode. It is the execution layer you
move to after the plain prescription is already correct.

- one prescription call like upstream `prop_run(...)`:
  - stay on `prop_run(...)`
- same prescription shape run repeatedly:
  - use `prepare_prescription(...)`
- same prescription run repeatedly or in parallel with different `PASSVALUE`s:
  - use `prepare_prescription_batch(...)`
- same prescription plus cached assets, slot-local state, or a named execution
  object:
  - use `prepare_model(...)`
- wavelength sweep:
  - build a vector of prepared runs and call `prop_run_multi(runs)`

Concrete translation:

```julia
psf, sampling = prop_run(my_prescription, 0.55, 256; PASSVALUE=passvalue)
```

```julia
prepared = prepare_prescription(my_prescription, 0.55, 256)
psf, sampling = prop_run(prepared; PASSVALUE=passvalue)
```

```julia
runs = [
    prepare_prescription(my_prescription, 0.50, 256; precision=Float32),
    prepare_prescription(my_prescription, 0.55, 256; precision=Float32),
    prepare_prescription(my_prescription, 0.60, 256; precision=Float32),
]

stack, samplings = prop_run_multi(runs)
```

For GPU-oriented throughput work, that last form is the intended public sweep
surface.

## Known Semantics Adopted In Julia
- `prop_resamplemap`: independent `xshift`/`yshift` handling
- `prop_end`: integer-safe centered extraction semantics
- `prop_state`: full wavefront state restoration in-place
- `prop_psd_errormap`: file reuse/write semantics and parity quirk handling are
  explicitly implemented
- `prop_rotate`: defaults to MATLAB-style linear interpolation; request cubic
  explicitly with `METH="cubic"` when needed
- `prop_magnify`: defaults to the damped-sinc `prop_szoom` path; request
  `QUICK=true` for cubic interpolation
- `prop_pixellate`: public API follows upstream PROPER PSF-integration
  semantics (`image, sampling_in, sampling_out, n_out`)
- `prop_writemap`: exported `XC_PIX` / `YC_PIX` FITS metadata follows the
  PROPER/MATLAB `floor(n/2)+1` convention

## Common Translation Notes
- wavelength inputs to `prop_run` stay in microns at the public boundary
- internal wavefront state uses meters
- mutating `...!` variants are available for callers that want explicit output
  ownership
- `prop_run_multi` is thread-based by default and preserves input ordering
- Julia does not expose Python-style runtime compatibility flags; accepted
  behavior choices are documented directly in [compatibility decisions](compat_decisions.md)

## Current Reference Workload
- The WFIRST Phase B reference port is used as a broad parity and benchmarking
  workload.
- It is a correctness and comparability target for the PROPER port, not a
  separate optimized execution path.

## Parity Interpretation
- Parity acceptance is evaluated with combined relative and absolute metrics.
- Deep-null/high-contrast outputs use denominator-floored relative metrics plus
  absolute thresholds.
- Threshold policy lives in:
  - [parity_thresholds.md](parity_thresholds.md)
  - `test/parity/thresholds/example_metrics_thresholds.json`
