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
- Read the [latency benchmarking guide](LATENCY_BENCHMARKING.md) before
  interpreting prepared-call tail latency or configuring percentile sample
  counts
- Run `./scripts/benchmark_latency.sh cpu`, `amdgpu`, or `cuda` for a
  correctness-gated service-time distribution and raw HdrHistogram log
- Run `./scripts/benchmark_cpu_gpu.sh` for the Julia CPU/GPU comparison lane
- Run `./scripts/benchmark_thread_topology.sh` for a correctness-gated CPU
  Julia/FFTW topology matrix; set `PROPER_THREAD_TOPOLOGY_WORKLOAD=batch` for
  independent prepared-run throughput
- Ordinary benchmark-only Julia dependencies live in `bench/Project.toml`; the
  exact-revision-pinned HdrHistogram overlay and version-specific manifests live
  in `bench/latency/`
- The benchmark scripts automatically develop the local checkout when needed
  and compose those environments for latency runs
- Optional GPU benchmark dependencies are declared separately in `bench/cuda/`
  and `bench/amdgpu/`; instantiate the project for the device available on the
  host
- The resulting summary is written to `bench/reports/julia_cpu_gpu_summary.md`
- The lane also writes:
  - `bench/reports/julia_cpu_gpu_batch_throughput.csv`
  - `bench/reports/julia_cpu_gpu_precision_split.csv`
- The summary includes a `Synthetic Core Propagation Tail` section for the
  shared propagation benchmark harness
- The summary also includes:
  - `Prepared Batch Throughput`
  - `Precision Split`
- Thread-topology reports default to
  `bench/reports/julia_thread_topology_cpu_core.json` and
  `bench/reports/julia_thread_topology_cpu_batch.json`; every candidate is
  checked against a reference output before its timing is accepted
- Run `./scripts/profile_core_cpu_gpu.sh` to capture backend-specific text
  profiles for that shared core workload
- The driver uses any available Julia GPU backends through the checked-in
  `bench/cuda/` and `bench/amdgpu/` environments
- BenchmarkTools-based rows in that lane are warmed steady-state results, not
  cold-start / TTFx measurements
- The scheduled prepared-latency lane collects artifacts without a percentile
  regression gate because shared CI hardware is not a stable latency runner
- Timing reports are hardware snapshots, not portable expectations. The current
  local machine validates AMDGPU; CUDA is availability-gated because no CUDA
  device is present here.
- If a GPU backend or supported device is unavailable, the summary records that
  backend as skipped rather than failing the whole benchmark lane
- Active GPU implementation tracking lives in
  [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)
- Python-baseline parity and Python-vs-Julia benchmarks use the accepted patched
  PROPER 3.3.4 source snapshot. Because SourceForge removed the historical
  archive, `./scripts/setup_python_baseline.sh` reproducibly reconstructs that
  exact snapshot from a checksum-pinned 3.3.5 archive. Run the script or set
  `PYPROPER_ROOT=/path/to/proper_v3.3.4_python` before running
  `./scripts/benchmark_all.sh`
- Python parity and WFIRST comparisons require the three upstream native C
  interpolation kernels; the harness compiles them in temporary storage and
  fails instead of silently substituting SciPy interpolation

## CI And Validation Workflows
- `.github/workflows/ci.yml` is the always-on regression workflow for pushes
  and pull requests
- `CI` runs unit tests across supported platforms, coverage/Codecov, and
  lightweight Python-baseline parity
- `CI` also executes all example smoke runners, builds this documentation with
  Documenter doctests, and checks local Markdown targets on Julia 1.10 and 1.12
- `.github/workflows/validation.yml` owns heavier validation surfaces:
  benchmark reports, prepared latency distributions, and the WFIRST Phase B
  Python/Julia parity matrix
- Benchmark/WFIRST validation runs on pushes to `main`, on a weekly schedule,
  and manually via `workflow_dispatch`; prepared latency runs on the schedule
  or manual dispatch
- Manual `Validation` runs can disable benchmark reports, prepared latency, or
  WFIRST parity and can limit WFIRST with the `wfirst_cases` input
- The verified Python PROPER reconstruction, upstream `proper-models` checkout,
  and WFIRST public data compatibility root are fetched into CI caches rather
  than committed to the repository

## Dependency Updates
- `.github/dependabot.yml` checks the package, test, documentation, example,
  coverage, benchmark, and GPU Julia environments every week
- Julia updates for the same dependency are grouped across environments so
  their compatibility bounds stay aligned; the local path dependency on
  `Proper` is excluded
- GitHub Actions updates are grouped into a separate weekly pull request
- Dependabot pull requests must pass the normal CI and validation contracts;
  the version-specific latency manifests remain part of that review

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
- [Parity closure report](PARITY_CLOSURE.md) (dated historical snapshot)
- [Semantic reconciliation report](SEMANTIC_RECONCILIATION.md) (dated
  historical snapshot)
- [WFIRST Phase B config matrix](WFIRST_PHASEB_CONFIG_MATRIX.md)

## Contributor Workflow
- [Porting checklist](PORTING_CHECKLIST.md)

## Maintainer Notes And Historical Records
These documents are useful when working on backend tuning, performance
investigation, or understanding how the port evolved. They are not the primary
entry point for new users.

- [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)
- [Archived maintainer notes](https://github.com/DarrylGamroth/Proper.jl/tree/main/docs/archive)
