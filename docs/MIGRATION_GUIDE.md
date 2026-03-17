# Migration Guide: Python/MATLAB PROPER to Proper.jl

## Scope
This guide summarizes compatibility expectations for users migrating prescriptions from Python PROPER 3.3.4 or MATLAB 3.3.1 to `Proper.jl`.

## Stable Expectations
- Public routine naming remains familiar (`prop_*`).
- One-to-one filename mapping is preserved for traceability.
- Python 3.3.4 remains the executable parity baseline.
- Default plotting in examples uses `Plots.jl`.
- FITS IO uses `FITSIO.jl`.
- The public execution contract is still tuple-based:
  - `prop_run(...) -> (psf, sampling)`
  - `prop_run_multi(...) -> (stack, samplings)`

## Compatibility Baseline
- Parity-first behavior against the patched Python 3.3.4 executable baseline.
- MATLAB/manual remain semantic references for investigating suspected translation defects.

## Practical Migration Shape
- Most prescriptions can keep familiar public calls:
  - `prop_begin`
  - `prop_lens`
  - `prop_propagate`
  - `prop_end`
- Compatibility keywords remain accepted in uppercase form.
- Lowercase Julia keyword aliases are also accepted at the public boundary.
- Julia internals are not a line-by-line transliteration:
  - mutable globals are replaced by explicit contexts and workspaces
  - prepared execution objects are available for repeated runs

### Recommended Julia Migration Path
1. Port the prescription with the familiar `prop_*` calls first.
2. Confirm parity against the Python baseline.
3. Introduce prepared execution only when repeated runs need workspace or asset reuse.

For a fuller prescription-porting workflow, see
`docs/PRESCRIPTION_AUTHORING_GUIDE.md`.

### Prepared Execution
- `prepare_prescription(...)` caches normalized execution state for a single prescription.
- `prepare_prescription_batch(...)` adds reusable per-slot contexts for repeated or parallel runs.
- `prepare_model(...)` adds optional prepared assets on top of that execution state.
- These are performance and ownership tools, not a separate compatibility mode.
- See `docs/PREPARED_EXECUTION_GUIDE.md` for concrete usage patterns.

## Known Semantics Adopted In Julia
- `prop_resamplemap`: independent `xshift`/`yshift` handling (aligned with MATLAB/manual intent).
- `prop_end`: integer-safe centered extraction semantics.
- `prop_state`: full wavefront state restoration in-place.
- `prop_psd_errormap`: file reuse/write semantics and parity quirk handling are explicitly implemented.
- `prop_rotate`: defaults to MATLAB-style linear interpolation; request cubic explicitly with `METH="cubic"` when needed.
- `prop_magnify`: defaults to the damped-sinc `prop_szoom` path; request `QUICK=true` for cubic interpolation.
- `prop_pixellate`: public API follows upstream PROPER PSF-integration semantics (`image, sampling_in, sampling_out, n_out`) rather than a simple integer box-downsample helper.
- `prop_writemap`: exported `XC_PIX` / `YC_PIX` FITS metadata now follows the PROPER/MATLAB `floor(n/2)+1` convention.

## Common Translation Notes
- Wavelength inputs to `prop_run` stay in microns at the public boundary.
- Internal wavefront state uses meters.
- Mutating `...!` variants are available for callers that want explicit output
  ownership.
- `prop_run_multi` in Julia is thread-based by default and preserves input
  ordering.
- Julia does not expose Python-style runtime compatibility flags; accepted
  behavioral choices are documented directly in `docs/compat_decisions.md`.

## Current Reference Workload
- The WFIRST Phase B reference port is used as a broad parity and benchmarking
  workload.
- It is intended as a correctness and comparability target for the PROPER port,
  not as a separate optimized execution path.

## Parity Interpretation
- Parity acceptance is evaluated with combined relative + absolute metrics.
- Deep-null/high-contrast outputs use denominator-floored relative metrics plus absolute thresholds.
- Threshold policy:
  - `docs/parity_thresholds.md`
  - `test/parity/thresholds/example_metrics_thresholds.json`
