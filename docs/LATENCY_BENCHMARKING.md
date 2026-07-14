# Latency Benchmarking

The prepared-latency lane answers a narrow question: after compilation,
planning, and allocation have completed, what distribution of service times does
one preallocated `PreparedRun` produce on this machine?

It complements the `BenchmarkTools` microbenchmarks. Those remain the right tool
for typical steady-state runtime and allocation comparisons; the latency lane is
for supported p50/p90/p99/p99.9 observations, raw histogram artifacts, and
run-to-run variation.

## Correctness Gate

Correctness is a prerequisite for measurement. Before recording any latency,
the driver:

1. computes an allocating CPU prepared-execution oracle;
2. executes the target preallocated CPU or GPU case twice;
3. requires the returned array to be the caller-owned output on both calls;
4. compares the first result with the CPU oracle and the second result with the
   first; and
5. verifies output sampling.

The driver aborts before timing when any check fails. GPU-to-host transfers used
for these comparisons are outside the timed region. Float64 is the default;
backend- and precision-specific tolerances are recorded in the JSON result and
can be overridden explicitly for an investigation.

## Measurement Contract

Each sample begins immediately before `prop_run(prepared_run)`.

- On CPU, the sample ends when `prop_run` returns and the caller-owned output is
  ready.
- On CUDA and AMDGPU, the sample ends after an explicit backend synchronization,
  so asynchronous launch time is not mistaken for completed service time.
- Preparation, FFT planning, compilation, correctness checks, report building,
  cooldown, and histogram export are excluded.
- Histogram recording occurs after the sample clock stops. Observer overhead is
  measured separately and reported, never subtracted.
- Runs are closed-loop with one outstanding operation. A new operation starts
  only after the preceding operation and its histogram update complete.
- The harness applies no coordinated-omission correction because it does not
  claim to model an external arrival schedule.

This is a service-time distribution, not a FAST-SCC control-loop or detector
arrival model. Fixed-rate arrivals, bursts, deadlines, queueing, dropped frames,
and overload behavior belong in the application or `AdaptiveOpticsSim.jl`
benchmark that owns those policies.

The histogram is configured for at least 1 ns through 60 s with three
significant figures by default. Values outside its representable range fail
instead of being silently clamped. Garbage-collector counters and achieved
completion rate are recorded for every run.

## Running The Benchmarks

From the repository root:

```sh
./scripts/benchmark_latency.sh cpu
./scripts/benchmark_latency.sh amdgpu
./scripts/benchmark_latency.sh cuda
./scripts/benchmark_latency.sh all
```

Unavailable GPU backends produce an explicit skipped JSON report. They do not
cause the CPU or another available backend to fail.

A short functional smoke run is useful before collecting a distribution:

```sh
PROPER_LATENCY_GRID_N=64 \
PROPER_LATENCY_SAMPLES=5 \
PROPER_LATENCY_REPETITIONS=1 \
PROPER_LATENCY_WARMUP=1 \
./scripts/benchmark_latency.sh amdgpu
```

For a p99 distribution with at least 100 observations in each run's upper one
percent, use at least 10,000 samples per run:

```sh
PROPER_LATENCY_GRID_N=512 \
PROPER_LATENCY_SAMPLES=10000 \
PROPER_LATENCY_REPETITIONS=3 \
PROPER_LATENCY_WARMUP=5 \
PROPER_LATENCY_COOLDOWN_SECONDS=1 \
./scripts/benchmark_latency.sh cpu
```

The main controls are:

| Variable | Default | Meaning |
| --- | ---: | --- |
| `PROPER_LATENCY_SAMPLES` | `1000` | Samples in each repeated run |
| `PROPER_LATENCY_REPETITIONS` | `3` | Independently reported histograms |
| `PROPER_LATENCY_WARMUP` | `3` | Untimed prepared calls before measurement |
| `PROPER_LATENCY_COOLDOWN_SECONDS` | `0` | Untimed pause between runs |
| `PROPER_LATENCY_GRID_N` | `512` | Square propagation-grid size |
| `PROPER_LATENCY_PRECISION` | `float64` | `float64`/`fp64` or `float32`/`fp32` |
| `PROPER_LATENCY_FFTW_THREADS` | `1` | FFTW threads, set before planning |
| `PROPER_LATENCY_BLAS_THREADS` | `1` | BLAS threads |
| `PROPER_LATENCY_MIN_TAIL_OBSERVATIONS` | `100` | Required observations beyond a reported tail percentile |
| `PROPER_LATENCY_LOWEST_NS` | `1` | Lowest discernible histogram value |
| `PROPER_LATENCY_HIGHEST_NS` | `60000000000` | Highest trackable histogram value |
| `PROPER_LATENCY_SIGNIFICANT_FIGURES` | `3` | Histogram precision, from 0 through 5 |
| `PROPER_LATENCY_ATOL` | backend-specific | Correctness absolute tolerance |
| `PROPER_LATENCY_RTOL` | backend-specific | Correctness relative tolerance |
| `PROPER_LATENCY_REPORT` | backend-specific path | JSON output override |
| `PROPER_LATENCY_HISTOGRAM` | backend-specific path | HdrHistogram log output override |

Thread settings are part of the experiment, not package defaults. Create a new
process and new prepared context after changing FFTW threads. Use the thread
topology benchmark first when selecting Julia/FFTW topology for a large CPU.

## Percentile Support

The report emits a percentile only when the histogram contains enough samples
to support it. With the default requirement of 100 observations beyond the
percentile, each individual histogram needs at least:

| Percentile | Minimum samples |
| ---: | ---: |
| p50 | 2 |
| p90 | 1,000 |
| p99 | 10,000 |
| p99.9 | 100,000 |

Unsupported requests appear in `unsupported_percentiles`; they are not
presented as estimates. Both per-run and aggregate histograms apply this rule
to their own sample counts. Do not infer p99.9 from a 10- or 20-sample
microbenchmark.

## Artifacts And Reproducibility

Each successful backend writes:

- `bench/reports/julia_prepared_latency_<backend>_<precision>.json`, containing
  the contract, environment, correctness evidence, GC counters, observer
  overhead, per-run histograms, aggregate histogram, and completion rate; and
- `bench/reports/julia_prepared_latency_<backend>_<precision>.hlog`, containing
  the tagged raw HdrHistogram interval logs.

Generated reports are intentionally ignored by Git and uploaded as CI
artifacts. The benchmark-only `HdrHistogram.jl` dependency is isolated under
`bench/latency/` and pinned to an exact source revision in Julia-version-specific
manifests. It is not a package runtime dependency.

Test the deterministic harness without running an optical workload:

```sh
julia +1.10 --startup-file=no --project=bench/latency \
  bench/latency/test/runtests.jl
julia +1.12 --startup-file=no --project=bench/latency \
  bench/latency/test/runtests.jl
```

## Validation Matrix

The matrix stays intentionally smaller than the Cartesian product of every
backend, precision, grid, thread count, and Julia version:

| Area | Permutations | Evidence | Status | Notes |
| --- | --- | --- | --- | --- |
| Histogram contract | Julia 1.10, 1.12 | deterministic CI tests for timing boundaries, range failure, percentile support, aggregation, and log round-trip | Covered | independent of optical hardware |
| CPU prepared execution | Float64 | scheduled correctness-gated 512-grid collection plus local smoke | Covered | shared CI artifacts are not regression gates |
| CPU prepared execution | Float32 | local correctness smoke | Gap | add automated optical integration coverage when Float32 latency becomes a maintained deployment contract |
| AMDGPU prepared execution | Float64, Float32 | local CPU-oracle and repeatability smoke on the available AMD GPU | Gap | requires a maintained AMDGPU runner for durable CI coverage |
| CUDA prepared execution | Float64, Float32 | unavailable-backend report on this machine | Gap | numerical coverage requires CUDA hardware; never infer it from AMDGPU |
| Thread topology | Julia/FFTW/BLAS choices | separate correctness-gated topology benchmark | Covered | do not multiply every latency row by every topology |

Update this table when a backend runner or maintained precision contract is
added. A local pass is useful evidence, but it is not a substitute for durable
CI coverage on that device class.

The scheduled shared-host CI lane collects artifacts and enforces configuration
and numerical correctness, but deliberately has no percentile regression gate.
Machine-level latency gates require a stable dedicated runner, controlled power
and affinity, repeated baselines, and a documented tolerance before they are
credible.
