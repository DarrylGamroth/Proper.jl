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

| MATLAB / manual mental model | Python PROPER usage | Julia usage |
| --- | --- | --- |
| `prop_run(prescription, λ, n, optval)` | `prop_run(prescription, λ, n, PASSVALUE=...)` | `prop_run(prescription, λ, n; PASSVALUE=...)` |
| same prescription called repeatedly | same prescription called repeatedly | `prepare_prescription(...)` |
| repeated runs with varying optional inputs | repeated runs with varying `PASSVALUE`s | `prepare_prescription_batch(...)` |
| one configured reusable model object | ad hoc application wrapper around one prescription | `prepare_model(...)` |
| wavelength sweep | wavelength sweep loop | `prop_run_multi(runs)` where `runs` is a vector of prepared runs |
| uppercase compatibility keywords | uppercase compatibility keywords | uppercase and lowercase keywords both accepted |

The practical rule is:
- start exactly where upstream usage starts
- only introduce prepared execution when repeated work or cached state is the
  actual problem

## Prescription Signature Translation

The most common upstream prescription shapes map to Julia like this.

| MATLAB | Python | Julia |
| --- | --- | --- |
| `function [psf, sampling] = p(lam, n)` | `def p(lam, n):` | `function p(λm, n)` |
| `function [psf, sampling] = p(lam, n, optval)` | `def p(lam, n, PASSVALUE=None):` | `function p(λm, n; radius=0.5, use_dm=false)` |
| positional optional input is common | positional third argument is common | prefer explicit Julia keywords; `PASSVALUE` remains an adapter for upstream-style calls |

`prop_run` normalizes map-like `PASSVALUE=...` values into ordinary Julia
keywords before calling the prescription. Non-map values still use the legacy
positional/`PASSVALUE` compatibility path.

## Side-By-Side Example

| MATLAB | Python | Julia |
| --- | --- | --- |
| ```matlab\nfunction [psf, sampling] = simple_prescription(lam, n, optval)\n    wf = prop_begin(1.0, lam, n);\n    wf = prop_circular_aperture(wf, 0.5);\n    wf = prop_define_entrance(wf);\n    wf = prop_lens(wf, 10.0);\n    wf = prop_propagate(wf, 10.0);\n    [psf, sampling] = prop_end(wf);\nend\n``` | ```python\ndef simple_prescription(lam, n, PASSVALUE=None):\n    wf = proper.prop_begin(1.0, lam, n)\n    proper.prop_circular_aperture(wf, 0.5)\n    proper.prop_define_entrance(wf)\n    proper.prop_lens(wf, 10.0)\n    proper.prop_propagate(wf, 10.0)\n    return proper.prop_end(wf)\n``` | ```julia\nusing Proper\n\nfunction simple_prescription(λm, n)\n    wf = prop_begin(1.0, λm, n)\n    prop_circular_aperture(wf, 0.5)\n    prop_define_entrance(wf)\n    prop_lens(wf, 10.0)\n    prop_propagate(wf, 10.0)\n    return prop_end(wf)\nend\n\npsf, sampling = prop_run(simple_prescription, 0.55, 256)\n``` |

What changes:
- Julia uses normal function definitions instead of Python module-level state.
- map-like `PASSVALUE` values are compatibility input and become Julia keywords.
- The optical sequence stays almost identical.

What does not change:
- wavelength at the public boundary is still in microns
- the prescription still returns `prop_end(wf)` or `(psf, sampling)`
- the familiar `prop_*` surface remains the main way to write the prescription

## DM / FITS / Map-Driven Prescription Pattern

This is the common nontrivial migration shape many upstream users care about.

| MATLAB / upstream-style idea | Python-style idea | Julia |
| --- | --- | --- |
| read FITS map, apply error map, optionally apply DM, propagate, end | read FITS map, apply error map, optionally apply DM, propagate, end | `migration_dm_fits_prescription` in [`examples/migration_dm_fits.jl`](../examples/migration_dm_fits.jl) |

The Julia version is:

```julia
using Proper

function dm_map_prescription(λm, n; errormap_path=nothing, dm_map=nothing)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)

    if errormap_path !== nothing
        prop_errormap(wf, errormap_path, WAVEFRONT=true, SAMPLING=prop_get_sampling(wf))
    end

    if dm_map !== nothing
        prop_dm(wf, dm_map)
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
    errormap_path="phase_map.fits",
    dm_map=zeros(256, 256),
)
```

The same prescription still accepts upstream-style map input through the
adapter:

```julia
psf, sampling = prop_run(
    dm_map_prescription,
    0.55,
    256;
    PASSVALUE=Dict("errormap_path" => "phase_map.fits", "dm_map" => zeros(256, 256)),
)
```

Important behavior:
- FITS files are read on the host, then promoted to the wavefront backend where
  feasible.
- if the FITS header does not carry pixel scale metadata, pass
  `SAMPLING=prop_get_sampling(wf)` or the appropriate map sampling explicitly
- `prop_dm(wf, dm_map)` expects a map already sampled on the wavefront grid.
- Use `prop_dm` with the full upstream-compatible actuator-space form only when
  you are intentionally modeling actuator geometry and influence functions.

## Keyword Translation Rules

- Uppercase compatibility keywords are accepted.
- Lowercase Julia aliases are also accepted at the public boundary.
- If you are writing new Julia-native code, prefer ordinary Julia keywords.

| MATLAB / manual style | Python style | Julia |
| --- | --- | --- |
| `prop_rotate(img, ang, 'METH', 'cubic')` | `proper.prop_rotate(img, ang, METH='cubic')` | `prop_rotate(img, ang; METH="cubic")` or `prop_rotate(img, ang; meth="cubic")` |
| `prop_magnify(img, mag, ..., 'QUICK')` | `proper.prop_magnify(img, mag, QUICK=True)` | `prop_magnify(img, mag; QUICK=true)` or `prop_magnify(img, mag; quick=true)` |
| `prop_errormap(..., 'WAVEFRONT')` | `proper.prop_errormap(..., WAVEFRONT=True)` | `prop_errormap(...; WAVEFRONT=true)` or `prop_errormap(...; wavefront=true)` |

## `PASSVALUE` Guidance

Use `PASSVALUE` when:
- porting an upstream prescription
- matching an existing parity harness
- carrying a heterogeneous set of optional inputs through `prop_run`

Prefer explicit Julia keywords or prepared assets when:
- the model is now Julia-native
- the same assets are reused across runs
- the call surface is stable enough to type explicitly

For new Julia-native prescriptions, prefer:

```julia
function run_coronagraph(λm, n; use_errors=false, occulter=:gaussian)
    # ...
end
```

over a dictionary-only signature. A compatibility wrapper may still translate
`PASSVALUE=Dict("use_errors" => true, "occulter_type" => "GAUSSIAN")` into
those keywords at the `prop_run` boundary.

Use symbols for user-facing selectors (`:gaussian`, `:solid`) and normalize to
typed singleton values internally when dispatch clarifies the implementation.

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

| MATLAB / Python intent | Julia |
| --- | --- |
| one run | `psf, sampling = prop_run(my_prescription, 0.55, 256; PASSVALUE=passvalue)` |
| one fixed prescription reused | `prepared = prepare_prescription(my_prescription, 0.55, 256)` then `prop_run(prepared; PASSVALUE=passvalue)` |
| wavelength sweep | `runs = [prepare_prescription(...), ...]` then `stack, samplings = prop_run_multi(runs)` |

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
