# Documentation Index

This directory contains both current user-facing documentation and historical
planning/audit records from the port effort.

## Start Here
- [API contract](api_contract.md)
- [Migration guide](MIGRATION_GUIDE.md)
- [Prescription authoring guide](PRESCRIPTION_AUTHORING_GUIDE.md)
- [Prepared execution guide](PREPARED_EXECUTION_GUIDE.md)
- [Runnable API examples](API_EXAMPLES.md)

## Benchmarking
- Run `./scripts/benchmark_cpu_gpu.sh` for the Julia CPU/GPU comparison lane
- The resulting summary is written to `bench/reports/julia_cpu_gpu_summary.md`
- The lane also writes:
  - `bench/reports/julia_cpu_gpu_batch_throughput.csv`
  - `bench/reports/julia_cpu_gpu_precision_split.csv`
- The summary includes a `Synthetic Core Propagation Tail` section for the
  shared propagation benchmark harness
- The summary also includes:
  - `Prepared Batch Throughput`
  - `Precision Split`
- Run `./scripts/profile_core_cpu_gpu.sh` to capture backend-specific text
  profiles for that shared core workload
- The driver uses any available Julia GPU backends (`CUDA.jl` and/or
  `AMDGPU.jl`)
- BenchmarkTools-based rows in that lane are warmed steady-state results, not
  cold-start / TTFx measurements
- If a GPU backend or supported device is unavailable, the summary records that
  backend as skipped rather than failing the whole benchmark lane
- Active GPU implementation tracking lives in
  [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)

## Core Contracts
- [Numerics contract](numerics_contract.md)
- [Backend traits contract](backend_traits.md)
- [Parity harness contract](parity_harness_contract.md)
- [Parity thresholds](parity_thresholds.md)
- [Compatibility decisions](compat_decisions.md)

## Reference Validation And Coverage
- [Parity closure report](PHASE8_CLOSURE.md)
- [Semantic reconciliation report](PHASE9_RECONCILIATION.md)
- [WFIRST Phase B config matrix](WFIRST_PHASEB_CONFIG_MATRIX.md)

## Contributor Workflow
- [Porting checklist](PORTING_CHECKLIST.md)

## Working Notes And Historical Records
These documents are useful when working on backend tuning, performance
investigation, or understanding how the port evolved. They are not the primary
entry point for new users.

- [Implementation progress](IMPLEMENTATION_PROGRESS.md)
- [GPU implementation plan](GPU_IMPLEMENTATION_PLAN.md)
- [Runtime optimization plan](RUNTIME_OPTIMIZATION_PLAN.md)
- [Performance follow-up plan](PERFORMANCE_FOLLOWUP_PLAN.md)
- [CUDA optimization plan](CUDA_OPTIMIZATION_PLAN.md)
- [Julia implementation audit (historical)](JULIA_IMPLEMENTATION_AUDIT_2026-03-04.md)
