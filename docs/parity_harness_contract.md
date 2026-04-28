# Parity Harness Contract

## Purpose
Define how Julia outputs are compared against Python 3.3.4 baselines and how artifacts/provenance are managed.

## Decision Links
- `D-0001` executable parity baseline
- `D-0010` golden data workflow (accepted)
- `D-0011` RNG determinism policy (accepted)
- `D-0014` CI/support matrix (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Baseline workflow accepted
- [x] CI jobs wired

## 1. Baseline Source
- Executable baseline: Python PROPER 3.3.4 from SourceForge:
  `https://sourceforge.net/projects/proper-library/files/proper_v3.3.4_python.zip`
- Local default path: `../proper_v3.3.4_python`
- Override path: `PYPROPER_ROOT=/path/to/proper_v3.3.4_python`
- Semantic tie-breakers: MATLAB 3.3.1 + PROPER manual

## 2. Artifact Layout
Recommended structure:
- `test/parity/cases/` case definitions (params, mode, backend)
- `test/parity/baseline/python334/` generated baseline arrays + metadata
- `test/parity/reports/` comparison summaries

## 3. Provenance Metadata Schema
Each baseline artifact should capture:
- python source snapshot identifier (commit/hash or equivalent)
- script/example name
- wavelength/grid/passvalue/config
- RNG seed info (if applicable)
- generation timestamp
- baseline identifier
- backend label and numeric precision

## 4. Metrics And Thresholds
Default gating approach:
- combine relative and absolute thresholds
- use denominator-floored relative metrics for near-zero/high-contrast outputs
- allow documented case-specific overrides

Additional metrics:
- max absolute error with case-specific bounds
- wrapped phase MAE bounds when phase comparisons apply

Executable threshold source:
- `test/parity/thresholds/example_metrics_thresholds.json`

## 5. Case Matrix
Minimum coverage:
- all 23 examples
- output mode (`NOABS=true/false`)
- map pipeline cases (`psdtest`, errormap/resample/rotate)
- multi-run ordering and shape semantics
- compat-mode divergence cases

## 6. Update Workflow
- Trigger baseline regeneration when:
  - accepted compatibility decision changes behavior
  - numerics contract changes
  - baseline generation scripts change
- Require report diff review before accepting updated baselines.

## 7. CI Contract
Required jobs:
- CPU unit tests on Linux, macOS, and Windows for supported Julia versions
- Linux coverage job
- Linux parity and benchmark job using the cached SourceForge Python PROPER
  baseline
Scheduled/extended jobs:
- full CPU parity and benchmark artifacts from the Linux parity job
Optional jobs:
- self-hosted CUDA tests
- self-hosted AMDGPU tests

## 8. Contract Tests
- [ ] Baseline generation reproducibility test (seeded)
- [ ] Artifact metadata completeness test
- [x] Threshold enforcement tests (`compare_examples.jl` + threshold config)
- [ ] CI job pass/fail integration test
