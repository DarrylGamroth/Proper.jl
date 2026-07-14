# Benchmark Suites

This directory holds benchmark-suite helpers that are shared by focused
performance measurements. Suite helpers should stay here unless they become
general infrastructure suitable for `bench/common/`.

Benchmark entry points currently include:

- `bench/julia/steady_state/dm_projection.jl`
- `bench/julia/steady_state/zernike_synthesis.jl`
- `bench/julia/steady_state/zernike_fit.jl`
- `bench/julia/steady_state/batch_preallocated_fp64.jl`
- `bench/julia/steady_state/thread_topology_case.jl` (normally launched through
  `scripts/benchmark_thread_topology.sh`)
- `bench/julia/wfirst_phaseb/prepared_models.jl`

Generated reports belong under ignored `bench/reports/`.

All reports include Julia, FFTW, and BLAS thread counts. Set
`PROPER_FFTW_THREADS` before launching a benchmark when comparing FFTW thread
configurations, and keep outer Julia-thread and inner FFTW-thread experiments
separate so oversubscription is visible.

Run the fresh-process topology matrix with:

```sh
PROPER_THREAD_TOPOLOGY_WORKLOAD=core \
  PROPER_THREAD_TOPOLOGY_MATRIX="1:1 1:4 1:8 1:16" \
  PROPER_BENCH_GRID_N=512 PROPER_BENCH_SAMPLES=20 \
  scripts/benchmark_thread_topology.sh
```

Each entry is `JULIA_THREADS:FFTW_THREADS`. BLAS defaults to one thread for
this FFT-dominated workload and can be changed with
`PROPER_BENCH_BLAS_THREADS`. Every case runs in a fresh process and must match
a common serialized field reference within the configured numerical tolerance
before its timing is retained. Set `PROPER_THREAD_TOPOLOGY_WORKLOAD=batch` to
measure the caller-owned four-run `prop_run_multi!` workload and compare outer
Julia threading, single-run FFTW threading, and hybrid configurations. Combined
reports are written to `bench/reports/julia_thread_topology_cpu_core.json` or
`bench/reports/julia_thread_topology_cpu_batch.json` unless
`PROPER_BENCH_REPORT` is set. The batch defaults to four wavelengths; set
`PROPER_BENCH_BATCH_SIZE` to a representative sweep size when evaluating a
many-core host.
