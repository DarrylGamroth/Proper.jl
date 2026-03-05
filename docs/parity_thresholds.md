# Parity Thresholds

This document defines executable parity acceptance thresholds for Python 3.3.4 baseline comparisons.

## Decision Links
- `D-0001`: Python 3.3.4 executable baseline
- `D-0030`: physics-equivalent parity goal
- `D-0031`: idiomatic Julia implementation style

## Policy
- Use both relative and absolute criteria.
- For deep-null/high-contrast outputs, absolute error and denominator-floored relative error are authoritative.
- Sampling parity is always required.
- Case-specific overrides are allowed and must be documented in threshold config.

## Threshold Source
- Machine-readable thresholds live in:
  - `test/parity/thresholds/example_metrics_thresholds.json`
- Enforcement is implemented in:
  - `test/parity/compare_examples.jl`

## Notes
- Relative-only thresholds are insufficient for coronagraph null regions where baseline values approach zero.
- Threshold updates require:
  - updated report evidence in `test/parity/reports/`
  - explicit review rationale in commit message and/or decision log
