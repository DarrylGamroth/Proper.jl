# Documentation Index

This directory contains both current user-facing documentation and historical
planning/audit records from the port effort.

## Start Here
- If you already know upstream PROPER:
  1. [Migration guide](MIGRATION_GUIDE.md)
  2. [Runnable API examples](API_EXAMPLES.md)
  3. [Prepared execution guide](PREPARED_EXECUTION_GUIDE.md)
- If you are writing or restructuring prescriptions in Julia after the initial
  port:
  - [Prescription authoring guide](PRESCRIPTION_AUTHORING_GUIDE.md)
- If you are integrating with AO/HIL runtimes such as `AdaptiveOpticsSim.jl`:
  - [Adaptive optics integration guide](ADAPTIVE_OPTICS_INTEGRATION.md)
- For the precise stable public surface:
  - [API contract](api_contract.md)

## Benchmarking
- Run `./scripts/benchmark_cpu_gpu.sh` for the Julia CPU/GPU comparison lane
- Benchmark-only Julia dependencies live in `bench/Project.toml`; the benchmark
  scripts automatically develop the local checkout into that environment before
  running benchmark code
- Optional GPU benchmark packages should be added to that environment, for
  example `julia --project=bench -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add("CUDA")'`
  or the same command with `"AMDGPU"`
- The resulting summary is written to `bench/reports/julia_cpu_gpu_summary.md`
- The lane also writes:
  - `bench/reports/julia_cpu_gpu_batch_throughput.csv`
  - `bench/reports/julia_cpu_gpu_precision_split.csv`
- The summary includes a `Synthetic Core Propagation Tail` section for the
  shared propagation benchmark harness
- The summary also includes:
  - `Prepared Batch Throughput`
  - `Precision Split`
- A representative validated CUDA batch result is now surfaced in the root
  README to make the prepared batch + FP32 capability visible from the main
  package entry point
- Run `./scripts/profile_core_cpu_gpu.sh` to capture backend-specific text
  profiles for that shared core workload
- The driver uses any available Julia GPU backends; CUDA or AMDGPU lanes require
  `CUDA.jl` or `AMDGPU.jl` to be available in the `bench/` environment
- BenchmarkTools-based rows in that lane are warmed steady-state results, not
  cold-start / TTFx measurements
- If a GPU backend or supported device is unavailable, the summary records that
  backend as skipped rather than failing the whole benchmark lane
- Active GPU implementation tracking lives in
  [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)
- Python-baseline parity and Python-vs-Julia benchmarks use PROPER 3.3.4 from
  SourceForge; run `./scripts/setup_python_baseline.sh` or set
  `PYPROPER_ROOT=/path/to/proper_v3.3.4_python` before running
  `./scripts/benchmark_all.sh`

## CI And Validation Workflows
- `.github/workflows/ci.yml` is the always-on regression workflow for pushes
  and pull requests
- `CI` runs unit tests across supported platforms, coverage/Codecov, and
  lightweight Python-baseline parity
- `.github/workflows/validation.yml` owns heavier validation surfaces:
  benchmark reports and the WFIRST Phase B Python/Julia parity matrix
- `Validation` runs on pushes to `main`, on a weekly schedule, and manually via
  `workflow_dispatch`
- Manual `Validation` runs can disable benchmark reports or WFIRST parity and
  can limit WFIRST with the `wfirst_cases` input
- The external Python PROPER baseline, upstream `proper-models` checkout, and
  WFIRST public data compatibility root are fetched into CI caches rather than
  committed to the repository

## Coverage
- Run `./scripts/coverage_lcov.sh` to execute the full package test suite with
  Julia coverage enabled and write `lcov.info`
- CI uploads `lcov.info` to Codecov using GitHub OIDC (`id-token: write` and
  `use_oidc: true`) instead of a shared `CODECOV_TOKEN`

## Core Contracts
- [Numerics contract](numerics_contract.md)
- [Backend traits contract](backend_traits.md)
- [Parity harness contract](parity_harness_contract.md)
- [Parity thresholds](parity_thresholds.md)
- [Compatibility decisions](compat_decisions.md)

## Reference Validation And Coverage
- [Parity closure report](PARITY_CLOSURE.md)
- [Semantic reconciliation report](SEMANTIC_RECONCILIATION.md)
- [WFIRST Phase B config matrix](WFIRST_PHASEB_CONFIG_MATRIX.md)

## Contributor Workflow
- [Porting checklist](PORTING_CHECKLIST.md)

## Maintainer Notes And Historical Records
These documents are useful when working on backend tuning, performance
investigation, or understanding how the port evolved. They are not the primary
entry point for new users.

- [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)
- [Archived maintainer notes](archive/)
