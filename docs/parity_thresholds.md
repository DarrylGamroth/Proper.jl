# Parity Thresholds

This document defines executable parity acceptance thresholds for Python 3.3.4 baseline comparisons.

## Decision Links
- `D-0001`: Python 3.3.4 executable baseline
- `D-0030`: physics-equivalent parity goal
- `D-0031`: idiomatic Julia implementation style
- `D-0060`: opt-in carrier-phase parity

## Policy
- Use both relative and absolute criteria.
- For deep-null/high-contrast outputs, absolute error and denominator-floored relative error are authoritative.
- Sampling parity is always required.
- Case-specific overrides are allowed and must be documented in threshold config.

## Threshold Source
- Machine-readable thresholds live in:
  - `test/parity/thresholds/example_metrics_thresholds.json`
  - `test/parity/cases/simple_case.toml` for the full-array simple case
  - `test/parity/cases/wavefront_accessors_even.toml` for centered accessors
  - `test/parity/cases/carrier_phase.toml` for coherent carrier tracking
- Enforcement is implemented in:
  - `test/parity/compare.jl`
  - `test/parity/compare_examples.jl`

## Full-Array Simple Case
- Relative L2 error must be at most `1e-12`.
- Sampling relative error must be at most `1e-12`.
- Non-finite or missing metrics fail the gate.
- CI runs the comparison with four Julia threads from the `examples`
  environment so the parity-only JSON/TOML dependencies remain outside the
  runtime package.

## Coherent Carrier Case

- Quarter-wave enabled/disabled and half-wave differential-arm complex fields
  use a maximum absolute error of `1e-14`.
- The destructive-interference mean intensity uses an absolute error of
  `1e-28`; an absolute gate is required because the expected value is near
  machine zero.

## Notes
- Relative-only thresholds are insufficient for coronagraph null regions where baseline values approach zero.
- Threshold updates require:
  - updated report evidence in `test/parity/reports/`
  - explicit review rationale in commit message and/or decision log
