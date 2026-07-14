# Parity Harness Contract

## Purpose
Define how Julia outputs are compared against Python 3.3.4 baselines and how artifacts/provenance are managed.

## Decision Links
- `D-0001` executable parity baseline
- `D-0010` golden data workflow (accepted)
- `D-0011` RNG determinism policy (accepted)
- `D-0014` CI/support matrix (accepted)
- `D-0070` tiered example validation evidence (accepted)

## Status
- [x] Draft pre-filled from proposed defaults
- [x] Baseline workflow accepted
- [x] CI jobs wired

## 1. Baseline Source
- Executable baseline: Python PROPER 3.3.4 from SourceForge:
  `https://sourceforge.net/projects/proper-library/files/proper_v3.3.4_python.zip`
- Local default path: `../proper_v3.3.4_python`
- Override path: `PYPROPER_ROOT=/path/to/proper_v3.3.4_python`
- External WFIRST/Roman Phase B Python model baseline:
  `https://github.com/ajeldorado/proper-models.git`
- Local default path: `../proper-models`
- Override path: `WFIRST_MODELS_ROOT=/path/to/proper-models` or
  `WFIRST_MODELS_PYTHON_ROOT=/path/to/wfirst_cgi/models_phaseb/python`
- CI uses a code-only `proper-models` checkout because the WFIRST Phase B Python
  model code does not require the repository's Git LFS assets. Set
  `WFIRST_MODELS_LFS=1` for a full local checkout when needed.
- Semantic tie-breakers: MATLAB 3.3.1 + PROPER manual

## 2. Artifact Layout
Current structure:
- `test/parity/cases/` case definitions (params, mode, backend)
- `test/parity/baseline/python334/` generated baseline arrays + metadata
- `test/parity/reports/` comparison summaries

## 3. Provenance Metadata Schema
Each metadata document uses schema version 1 and captures:
- baseline name/version/source URL, local source root, and deterministic SHA-256
  of the baseline source tree
- generator entry point and generation timestamp
- Python executable/version plus NumPy, SciPy, and Astropy versions
- backend label and numeric precision
- RNG seed (or an explicit null value)
- per-case wavelength, grid size, and configuration
- SHA-256 for every generated artifact

`test/parity/validate_artifacts.jl` enforces the schema and artifact hashes.

## 4. Metrics And Thresholds
Default gating approach:
- combine relative and absolute thresholds
- use denominator-floored relative metrics for near-zero/high-contrast outputs
- allow documented case-specific overrides

Additional metrics:
- max absolute error with case-specific bounds
- wrapped phase MAE bounds when phase comparisons apply
- center values and four asymmetric pixel probes for example matrices, including
  real and imaginary components for complex fields; these catch transposes,
  reflections, and centering errors that aggregate norms can hide

Executable threshold source:
- `test/parity/thresholds/example_metrics_thresholds.json`

## 5. Case Matrix
Minimum coverage is tiered by what each upstream source actually exposes:
- filename, load, and expected-entry-point coverage for all 23 upstream example
  source files: one package marker plus 22 executable/helper/demo modules
- headless reduced-grid CPU execution for every user-facing runner
- numerical parity at meaningful prescription/runner output boundaries; the
  current executable matrix contains 16 cases, including broadband
  `testmulti1` and all three `testmulti2` ripple patterns
- helper and plotting wrappers covered transitively or by execution smoke when
  they do not expose a distinct numerical return
- output mode (`NOABS=true/false`)
- map pipeline cases (`psdtest`, errormap/resample/rotate)
- multi-run ordering and shape semantics
- accepted or corrected behavior cases recorded in `docs/compat_decisions.md`

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
- Linux Python-baseline parity job using the cached SourceForge Python PROPER
  baseline
Scheduled/extended jobs:
- benchmark artifact generation from the full benchmark suite
- WFIRST Phase B Python/Julia parity matrix using the cached SourceForge Python
  PROPER baseline, cached upstream `proper-models` checkout, and cached public
  Roman preflight compatibility data
- the extended validation workflow runs on pushes to `main`, weekly schedule,
  and manual `workflow_dispatch`; manual runs can target selected WFIRST cases
Optional jobs:
- self-hosted CUDA tests
- self-hosted AMDGPU tests

## 8. Contract Tests
- [x] Baseline generation reproducibility test (seeded and byte-for-byte via
  `scripts/check_parity_reproducibility.sh`)
- [x] Artifact metadata completeness and SHA-256 test
  (`test/parity/validate_artifacts.jl`)
- [x] Threshold enforcement tests (`compare_examples.jl` + threshold config)
- [x] CI integration for regeneration, reproducibility, metadata validation,
  and threshold enforcement (`.github/workflows/ci.yml`)
