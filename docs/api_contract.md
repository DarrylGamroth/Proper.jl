# API Contract

## Purpose
Define the user-facing API guarantees for `Proper.jl` so ports remain familiar to PROPER users while using idiomatic Julia internals.

## Decision Links
- `D-0001` baseline/parity policy
- `D-0004` package/module naming (accepted)
- `D-0005` public API compatibility surface (accepted)
- `D-0035` remove runtime compatibility modes (accepted)
- `D-0028` dynamic dispatch minimization (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Decisions accepted and reflected
- [ ] Tests aligned

## 1. Package And Module Names
- Canonical module name: `Proper`
- Compatibility alias strategy:
  - Keep PROPER-familiar `prop_*` function names.
  - Permit lightweight example-level aliasing (`const proper = Proper`) for user familiarity.
  - Add a separate wrapper module/package only if needed later.
- Export policy:
  - Export public `prop_*` routines intended for end users.
  - Export core public types (`WaveFront`, config/runtime types).

## 2. Public Entry Points
Stable entry points:
- `prop_run`
- `prop_run_multi`
- `prepare_prescription`
- `prepare_prescription_batch`
- `prepare_model`
- `prop_begin`
- `prop_end`
- `prop_propagate`
- `prop_lens`
- `prop_dm`
- `prop_errormap`
- `prop_psd_errormap`

Notes:
- Keep PROPER-style names for user familiarity.
- Internals may be refactored freely if behavior contract is preserved.
- Reusable runtime state may be supplied explicitly:
  - `prop_run(...; context=ctx)`
  - `prop_run(prepare_prescription(...))`
  - `prop_run(prepare_prescription_batch(...); slot=1)`
  - `prop_begin(...; context=ctx)` / `prop_begin(...; workspace=ws)`
  - `prop_wavefront(...; context=ctx)` / `prop_wavefront(...; workspace=ws)`
- Prepared parallel execution forks stored runtime state per pass:
  - `prop_run_multi(prepared::PreparedPrescription)` clones the prepared `RunContext` into independent workspaces before threaded execution.
  - This preserves backend/planning configuration without sharing mutable workspace scratch across threads.
  - `prop_run_multi(prepared_batch::PreparedBatch)` reuses a growable pool of those forked contexts across repeated calls.

## 3. Parity Baseline
- Behavior targets the patched Python 3.3.4 executable baseline used by the parity harness.
- MATLAB/manual references are used to evaluate suspected translation defects.

### 3.1 Context Constructor Contract
- Canonical constructor entry point (name subject to implementation details):
  - `RunContext(; backend=..., rng=..., workspace=..., fft_planning=...)`
- Prepared execution object:
  - `PreparedPrescription`
  - `PreparedBatch`
  - `PreparedModel`
  - `prepare_prescription(routine_name, lambda0_microns, gridsize; context=..., PASSVALUE=..., kwargs...)`
  - `prepare_prescription_batch(prepared_or_routine, ...; pool_size=...)`
  - `prepare_model(prepared_or_routine, ...; pool_size=..., assets=...)`
- `prop_*` public APIs should consume `RunContext` (or equivalent typed config) without compatibility mode flags.

## 4. Keyword Argument Contract
- Keyword style support:
  - [x] uppercase compatibility keywords (`NORM`, `NOABS`, ...)
  - [x] lowercase Julia aliases (`norm`, `noabs`, ...)

### 4.1 Canonicalization Rules
- Precedence if both forms provided:
  - if values are equivalent, accept
  - if values conflict, throw `ArgumentError`
- Unknown keyword handling:
  - throw `ArgumentError` with allowed-keyword hints
- Type coercions:
  - accept `Bool` natively
  - allow legacy `0/1` integers for boolean switches at compatibility boundary

## 5. Return Value Contract
- `prop_run(...) -> (psf, sampling)` where:
  - `psf`: `AbstractMatrix{<:Number}` (typically real intensity, complex when requested)
  - `sampling`: `Float64` (meters per pixel)
- `prop_end(...) -> (wavefront_or_intensity, sampling)` where:
  - `NOABS=false`: real intensity matrix
  - `NOABS=true`: complex field matrix
- `prop_run_multi(...) -> (stacked_psf, sampling_vector)` where:
  - `stacked_psf`: `AbstractArray{<:Number,3}`
  - `sampling_vector`: `Vector{Float64}`

## 6. Errors And Warnings
- Error type policy:
  - user input/keyword/config errors: `ArgumentError` / `DomainError`
  - missing backend implementation: `MethodError`
  - unsupported runtime path: `ErrorException` with actionable message
- Warning policy:
  - `@warn` for compatibility-relevant behavior changes

## 7. Non-Goals
- Reproducing Python module-global mutable state patterns.
- Enforcing Python object model or dynamic import behavior internally.
- Allowing runtime symbol/dict-based branching in hot inner loops.

## 8. Contract Tests
- [ ] API smoke tests for all stable entry points
- [ ] Keyword compatibility tests (uppercase/lowercase)
- [ ] Return shape/type tests
- [ ] Baseline parity behavior tests
