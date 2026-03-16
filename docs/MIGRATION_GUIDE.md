# Migration Guide: Python/MATLAB PROPER to Proper.jl

## Scope
This guide summarizes compatibility expectations for users migrating prescriptions from Python PROPER 3.3.4 or MATLAB 3.3.1 to `Proper.jl`.

## Stable Expectations
- Public routine naming remains familiar (`prop_*`).
- One-to-one filename mapping is preserved for traceability.
- Python 3.3.4 remains the executable parity baseline.
- Default plotting in examples uses `Plots.jl`.
- FITS IO uses `FITSIO.jl`.

## Compatibility Baseline
- Parity-first behavior against the patched Python 3.3.4 executable baseline.
- MATLAB/manual remain semantic references for investigating suspected translation defects.

## Known Semantics Adopted In Julia
- `prop_resamplemap`: independent `xshift`/`yshift` handling (aligned with MATLAB/manual intent).
- `prop_end`: integer-safe centered extraction semantics.
- `prop_state`: full wavefront state restoration in-place.
- `prop_psd_errormap`: file reuse/write semantics and parity quirk handling are explicitly implemented.
- `prop_rotate`: defaults to MATLAB-style linear interpolation; request cubic explicitly with `METH="cubic"` when needed.
- `prop_pixellate`: public API follows upstream PROPER PSF-integration semantics (`image, sampling_in, sampling_out, n_out`) rather than a simple integer box-downsample helper.
- `prop_writemap`: exported `XC_PIX` / `YC_PIX` FITS metadata now follows the PROPER/MATLAB `floor(n/2)+1` convention.

## Parity Interpretation
- Parity acceptance is evaluated with combined relative + absolute metrics.
- Deep-null/high-contrast outputs use denominator-floored relative metrics plus absolute thresholds.
- Threshold policy:
  - `docs/parity_thresholds.md`
  - `test/parity/thresholds/example_metrics_thresholds.json`
