# Numerics Contract

## Purpose
Freeze numerical conventions used by propagation, transforms, and coordinate systems to prevent silent drift.

## Decision Links
- `D-0001` baseline/parity policy
- `D-0035` remove runtime compatibility modes (accepted)
- `D-0007` numerical convention contract (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Conventions frozen
- [ ] Regression tests in place

## 1. Precision Policy
- Default floating precision: `Float64`
- Default complex field precision: `ComplexF64`
- Prepared execution precision:
  - `prepare_prescription`, `prepare_prescription_batch`, and `prepare_model`
    accept `precision=Float32` or `precision=Float64`
  - when explicit precision is requested, the prepared wavelength, workspace,
    and downstream field allocations follow that precision
  - if `context=...` is supplied alongside `precision=...`, the context
    workspace precision must match the requested precision
- GPU precision policy:
  - preserve input eltype/backend where practical
  - avoid implicit downcast unless explicitly requested

## 2. FFT Convention
- Backend primitive convention:
  - use `AbstractFFTs`/backend-default semantics (forward unscaled, inverse scaled by `1/N` product)
- PROPER routine convention:
  - apply explicit routine-level scaling to match the patched Python 3.3.4 baseline
- Net scaling expectation:
  - parity is validated at routine outputs, not raw backend FFT calls
- Backend equivalence tolerance:
  - must satisfy parity thresholds in plan (`rel L2` CPU/GPU)

## 3. Centering Convention
- `prop_shift_center` semantics:
  - forward shift uses `floor(size/2)` on each axis
  - inverse shift uses `ceil(size/2)` on each axis
- Even/odd behavior:
  - self-inverse guarantee applies for even dimensions
  - odd-dimension behavior is explicit and regression-tested
- Direction/inverse semantics:
  - explicit inverse variant may be provided for corrected/internal use

## 4. Coordinate And Pixel Center Convention
- Pixel-center definition:
  - coordinate zero is defined per-routine according to accepted upstream semantics
- Beam center index rule:
  - interpolation/resampling routines that follow MATLAB use `fix(n/2)+1` semantics
  - routines that intentionally preserve upstream Python/C executable behavior document that choice in `docs/compat_decisions.md`
- Coordinate orientation:
  - internal array math follows PROPER/Python indexing conventions
  - display orientation choices are handled in plotting layer

## 5. Units Contract
- Internal wavelength/sampling units: meters
- `prop_run` wavelength input: microns (converted to meters internally)
- Sampling outputs: meters per pixel
- Map unit conversions:
  - support meters, nanometers, microns with explicit switch/keyword handling

## 6. Randomness Contract
- Runtime default RNG behavior:
  - non-deterministic by default (familiar PROPER behavior)
- Test/parity seeded behavior:
  - explicit seed required
- Cross-backend parity workflow:
  - generate deterministic random phase maps on CPU and transfer as needed

## 7. Numerical Divergence Policy
- Preserve patched Python 3.3.4 baseline behavior unless a divergence is explicitly documented.
- Any divergence must include:
  - rationale tied to MATLAB/manual intent or correctness
  - parity evidence and updated thresholds where relevant
  - regression tests covering the chosen behavior

## 8. Contract Tests
- [ ] FFT normalization regression tests
- [ ] Centering tests (even and odd grid)
- [ ] Coordinate convention tests
- [ ] Unit conversion tests
- [ ] Seeded PSD/error-map reproducibility tests
