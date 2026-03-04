# Numerics Contract

## Purpose
Freeze numerical conventions used by propagation, transforms, and coordinate systems to prevent silent drift.

## Decision Links
- `D-0001` baseline/parity policy
- `D-0006` compat mode default (accepted)
- `D-0007` numerical convention contract (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Conventions frozen
- [ ] Regression tests in place

## 1. Precision Policy
- Default floating precision: `Float64`
- Default complex field precision: `ComplexF64`
- GPU precision policy:
  - preserve input eltype/backend where practical
  - avoid implicit downcast unless explicitly requested

## 2. FFT Convention
- Backend primitive convention:
  - use `AbstractFFTs`/backend-default semantics (forward unscaled, inverse scaled by `1/N` product)
- PROPER routine convention:
  - apply explicit routine-level scaling to match Python 3.3.4 in `:python334`
- Net scaling expectation:
  - parity is validated at routine outputs, not raw backend FFT calls
- Backend equivalence tolerance:
  - must satisfy parity thresholds in plan (`rel L2` CPU/GPU)

## 3. Centering Convention
- `prop_shift_center` semantics:
  - 2D circular shift by `floor(size/2)` on both axes (Python parity)
- Even/odd behavior:
  - self-inverse guarantee applies for even dimensions (primary PROPER path)
  - odd-dimension behavior is defined but non-primary
- Direction/inverse semantics:
  - explicit inverse variant may be provided for corrected/internal use

## 4. Coordinate And Pixel Center Convention
- Pixel-center definition:
  - coordinate zero at index corresponding to `n ÷ 2 + 1` in Julia (for even `n`)
- Beam center index rule:
  - follows Python expressions based on `(arange(n) - n//2)` equivalence
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

## 7. Compat-Mode Numerical Differences
Document expected differences:
- `:python334`:
  - preserve Python 3.3.4 behavior, including known quirks where intentionally retained
- `:corrected`:
  - apply documented fixes (e.g., map-shift semantics, extract indexing, state restore semantics, backend-toggle side effects)

## 8. Contract Tests
- [ ] FFT normalization regression tests
- [ ] Centering tests (even and odd grid)
- [ ] Coordinate convention tests
- [ ] Unit conversion tests
- [ ] Seeded PSD/error-map reproducibility tests
