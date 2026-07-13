# Numerics Contract

## Purpose
Freeze numerical conventions used by propagation, transforms, and coordinate systems to prevent silent drift.

## Decision Links
- `D-0001` baseline/parity policy
- `D-0035` remove runtime compatibility modes (accepted)
- `D-0007` numerical convention contract (accepted)
- `D-0057` centered public wavefront accessors (accepted)
- `D-0063` square propagation grid and independent transform oracle (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Conventions frozen
- [x] Regression tests in place for the centered-field and map-application
  contracts changed by `D-0057`

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
- Internal/public ordering:
  - `WaveFront.field` remains in raw FFT ordering for propagation and kernels
  - `prop_get_wavefront`, `prop_get_amplitude`, and `prop_get_phase` return
    centered backend-preserving arrays
  - matrix inputs to `prop_add_wavefront` are centered public data and are
    inverse-shifted before being added to the internal field
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

## 8. Correctness Coverage Matrix

| Area | Axis | Permutations | Evidence | Status | Notes |
| --- | --- | --- | --- | --- | --- |
| Public field ordering | Grid shape/parity | `7x7`, `8x8`, `7x8`, `8x7` | `test/test_public_helper_coverage.jl` | Covered | Asymmetric complex values catch wrong shifts, transposes, and aliases. |
| Centered field addition | Input form | scalar, centered matrix | `test/test_public_helper_coverage.jl` | Covered | Matrix extraction and reinsertion compose on odd and even axes. |
| Backend preservation | Backend/shape/layout | CPU; CUDA and AMDGPU on odd, even, mixed rectangular, and stepped `SubArray` inputs | `test/test_r2_trait_routing.jl` | Covered | GPU checks run when the corresponding CI runner is enabled, prohibit scalar indexing, and assert KA routing for stepped views. |
| FFT-scratch safety | Operation/backend | multiply and divide on CPU, CUDA, AMDGPU | `test/test_complex_map_scratch_alias.jl` | Covered | Regression reproduces a planned point-to-point propagation before applying an asymmetric complex map. |
| FFT planning equivalence | Policy/transform/direction | `ESTIMATE`, `MEASURE`; PTP/WTS/STW; positive and negative distance; live-scratch replan | `test/test_run_context_correctness.jl`, `test/test_r3_mutating_workspace.jl` | Covered | Numerical complex-field assertions prove planning cannot mutate live input. |
| Coherent carrier phase | Path/execution/backend | PTP/WTS/STW, split propagation, two-arm null; direct/prepared/batch/model/prepared-run/multi; CPU and optional CUDA/AMDGPU | `test/test_run_context_correctness.jl`, `test/test_r2_trait_routing.jl` | Covered | Default remains envelope-only; enabled mode asserts complex fields and half-wave destructive interference. |
| Propagation grid shape | Shape/entry point/state | square accepted; rectangular PTP/WTS/STW/`prop_propagate` rejected before mutation | `test/test_propagation_numerics.jl` | Covered | Physical propagation retains upstream's single-sampling square-grid contract; rectangular storage remains available to general matrix helpers. |
| FFT normalization oracle | Transform/direction/grid | direct-loop DFT vs PTP/WTS/STW, positive and negative distances, odd grid | `test/test_propagation_numerics.jl` | Covered | Oracle does not call FFTW; assertions cover complex field values, power conservation, PTP inverse roundtrip, and sampling updates. |
| Coordinate convention | Grid parity/order | even and odd spatial axes; centered and raw FFT frequency order | `test/test_propagation_numerics.jl` | Covered | Exact arrays catch half-pixel, shift-direction, and axis-order drift. |
| Unit conversion | API boundary/unit | run wavelength microns-to-meters, sampling radians/arcseconds, error maps in nm/microns | `test/test_propagation_numerics.jl` | Covered | Unit tests assert physical scale factors at public boundaries. |
| Prepared RNG ownership | RNG source/scheduling/planning | implicit context RNG, explicit override, serial slots, threaded batch, Estimate/Measure | `test/test_run_context_correctness.jl`; CI package tests use `-t4` | Covered | Identically seeded context trees produce exact ordered stacks; explicit RNG leaves context RNG untouched. |
| Executable upstream parity | Python 3.3.4 accessor/add and carrier paths | centered field, amplitude, phase, scalar add, matrix add; quarter-wave carrier and two-arm null | `test/parity/generate_python_baseline.py`, `test/parity/compare.jl` | Covered | Exact pixelwise threshold for field/add outputs plus complex-field and null-intensity thresholds for opt-in carrier tracking. |
| Threaded stack assembly | Packed output storage | `BitMatrix` outputs into `BitArray{3}` | `test/test_multi_run_scheduling.jl`; CI package tests use `-t4` | Covered | Repeated yields expose logical slices that share packed storage words. |

Additional contract areas retained for follow-up are tracked in the porting
plan as new compatibility surfaces are added.
