# Proper.jl

[![CI](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/ci.yml)
[![Validation](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/validation.yml/badge.svg)](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/validation.yml)
[![codecov](https://codecov.io/gh/DarrylGamroth/Proper.jl/graph/badge.svg)](https://codecov.io/gh/DarrylGamroth/Proper.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia port of John Krist's PROPER optical propagation library.

PROPER is a diffraction propagation toolkit for exploring wave-optics behavior
in optical systems. It is not a ray-tracing system and is not intended to
replace detailed optical design tools. `Proper.jl` keeps the familiar public
`prop_*` API while using Julia-native internals, explicit runtime state, and a
prepared execution layer for repeated runs.

## Start Here If You Already Know PROPER

If you are coming from PROPER for Python or MATLAB, the intended on-ramp is:

1. write the prescription in Julia with the familiar `prop_*` calls
2. run it with `prop_run(...)`
3. prefer explicit Julia keywords for stable model options
4. keep `PASSVALUE` only when you actually need upstream-style compatibility
5. move to prepared execution only after the plain prescription is correct

The one primary migration document is:
- [docs/MIGRATION_GUIDE.md](docs/MIGRATION_GUIDE.md)

The minimum supporting docs after that are:
- [docs/API_EXAMPLES.md](docs/API_EXAMPLES.md)
- [docs/PREPARED_EXECUTION_GUIDE.md](docs/PREPARED_EXECUTION_GUIDE.md)

You should not need to search the repository to figure out how to port a normal
upstream prescription.

## Project Status
- public `prop_*` routine naming preserved for migration familiarity
- Python PROPER 3.3.4 is the executable parity baseline
- MATLAB PROPER 3.3.1 is used as a semantic reference, especially for
  column-major correctness questions
- broad CPU parity coverage is in place, including WFIRST Phase B reference
  workloads
- prepared execution (`PreparedPrescription`, `PreparedBatch`,
  `PreparedModel`) is available for repeated and parallel runs

## Quick Migration Example

Upstream Python PROPER users usually start from a prescription like this:

```python
def simple_prescription(lam, n, PASSVALUE=None):
    wf = proper.prop_begin(1.0, lam, n)
    proper.prop_circular_aperture(wf, 0.5)
    proper.prop_define_entrance(wf)
    return proper.prop_end(wf)
```

The direct Julia shape is:

```julia
using Proper

function simple_prescription(λm, n; aperture_radius=0.5)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, aperture_radius)
    prop_define_entrance(wf)
    return prop_end(wf)
end

psf, sampling = prop_run(simple_prescription, 0.55, 128; aperture_radius=0.5)
```

That is the default migration path. Start there before introducing prepared
execution, explicit contexts, or GPU backends.

`PASSVALUE` is still accepted for upstream-style ports and parity workflows,
but Julia-native prescriptions should prefer ordinary keywords and symbols such
as `occulter=:gaussian` over string selector dictionaries.

## Running A Prescription

The familiar execution contract is preserved:

```julia
(psf, sampling) = prop_run(prescription, wavelength_microns, grid_size; kwargs...)
```

where:
- `prescription` is a function or globally resolvable name
- `wavelength_microns` is the wavelength in microns at the public boundary
- `grid_size` is the computational grid dimension

### Minimal Example
```julia
using Proper

function simple_prescription(λm, n)
    wf = prop_begin(1.0, λm, n)
    prop_circular_aperture(wf, 0.5)
    prop_define_entrance(wf)
    return prop_end(wf)
end

psf, sampling = prop_run(simple_prescription, 0.55, 128)
```

## Choose The Right Execution Surface

### Use `prop_run(...)` First
Use plain `prop_run` when:
- porting or validating a prescription
- comparing behavior against upstream Python or MATLAB references
- running a prescription occasionally

```julia
psf, sampling = prop_run(simple_prescription, 0.55, 128)
```

### Use `prepare_prescription(...)` When Runs Repeat
Use `PreparedPrescription` when the prescription shape stays fixed and you want
to reuse normalized execution state.

```julia
prepared = prepare_prescription(simple_prescription, 0.55, 128)
psf, sampling = prop_run(prepared)
```

You can request an explicit prepared execution precision when you want the
repeated execution surface to stay on `Float32` or `Float64`:

```julia
prepared = prepare_prescription(simple_prescription, 0.55f0, 128; precision=Float32)
psf, sampling = prop_run(prepared)
```

### Use `prepare_prescription_batch(...)` For Repeated Parallel Work
Use `PreparedBatch` when you want reusable per-slot contexts for repeated or
parallel execution.

```julia
batch = prepare_prescription_batch(simple_prescription, 0.55, 128; pool_size=2)
stack, samplings = prop_run_multi(batch; PASSVALUE=[nothing, nothing])
```

For wavelength sweeps or mixed prepared execution objects, pass a vector of
prepared runs directly:

```julia
runs = [
    prepare_prescription(simple_prescription, 0.50f0, 128; precision=Float32),
    prepare_prescription(simple_prescription, 0.55f0, 128; precision=Float32),
    prepare_prescription(simple_prescription, 0.60f0, 128; precision=Float32),
]

stack, samplings = prop_run_multi(runs)
```

This vector-of-prepared-runs form is the intended throughput surface for fixed
wavelength sweeps. It is the same API used by the batch-throughput benchmark
lane.

### Use `prepare_model(...)` When Assets Or Naming Matter
Use `PreparedModel` when you want a named execution object and optional cached
assets layered on top of prepared contexts.

For users coming from upstream PROPER, `PreparedModel` should be read as the
Julia reusable execution object for one configured prescription, not as a
different optical model type.

```julia
model = prepare_model(:simple_model, simple_prescription, 0.55, 128; pool_size=2)
psf, sampling = prop_run(model; slot=1)
```

For a fuller explanation of prepared execution, see
`docs/PREPARED_EXECUTION_GUIDE.md`.

## Compatibility Notes
- Python PROPER 3.3.4 is the executable parity baseline
- MATLAB PROPER 3.3.1 remains the semantic reference when investigating likely
  translation defects
- accepted behavior choices are recorded in `docs/compat_decisions.md`
- no runtime compatibility-mode flags are exposed; behavior decisions are
  documented and tested directly

Important current examples:
- `prop_rotate` defaults to linear interpolation; request cubic explicitly with
  `METH="cubic"` or `CUBIC=true`
- `prop_magnify` defaults to the damped-sinc `prop_szoom` path; `QUICK=true`
  selects cubic interpolation
- `prop_pixellate` follows the upstream PROPER detector-integration API:
  `prop_pixellate(image, sampling_in, sampling_out, n_out=0)`

## Examples And Reference Workloads
- `examples/` contains ported example prescriptions and smoke examples
- [`examples/migration_dm_fits.jl`](examples/migration_dm_fits.jl) is the
  concrete migration example for a `PASSVALUE` + FITS errormap + DM map
  workflow
- the WFIRST Phase B reference port is included as a broad correctness and
  benchmarking workload
- the WFIRST reference model lives under `reference_models/wfirst_phaseb_proper`

## Benchmark And Parity Summary
- Parity closure is complete against the patched Python 3.3.4 baseline:
  see [docs/PARITY_CLOSURE.md](docs/PARITY_CLOSURE.md)
- MATLAB/manual semantic reconciliation is complete for the identified hotspot
  semantics: see [docs/SEMANTIC_RECONCILIATION.md](docs/SEMANTIC_RECONCILIATION.md)
- The WFIRST Phase B reference port is used as a broad correctness and
  benchmarking workload, not as a separate optimized execution path
- Current benchmark policy separates steady-state runtime from Julia cold-start /
  TTFx (`D-0029`)

## CI And Validation Layout
Benchmarking and WFIRST validation stay in this repository for traceability,
but they are separated from the always-on package CI:

- `CI` runs package tests, coverage, Codecov upload, and lightweight
  Python-baseline parity on pushes and pull requests
- `Validation` runs benchmark reports and the heavy WFIRST Phase B parity matrix
  on pushes to `main`, weekly schedule, or manual `workflow_dispatch`
- manual `Validation` runs can disable either benchmark reports or WFIRST
  parity, and can override `wfirst_cases` for a targeted WFIRST subset
- external baselines are fetched into CI caches rather than vendored into the
  repository

This keeps normal pull request feedback focused on regressions while preserving
the full benchmark/WFIRST validation surface in the same commit history as the
implementation.

## Coverage
Run the local coverage lane with:

```bash
./scripts/coverage_lcov.sh
```

This runs the package tests with Julia coverage enabled and writes `lcov.info`.
CI uploads the same report to Codecov using GitHub OIDC, not a shared upload
secret.

## Validation Workflows
The full benchmark and WFIRST lanes are intentionally not required for every
pull request. To run them locally:

```bash
./scripts/setup_python_baseline.sh
./scripts/setup_wfirst_models_baseline.sh
./scripts/setup_parity_venv.sh

# Benchmark reports without WFIRST.
PYTHON_BIN=.venv-parity/bin/python BENCH_INCLUDE_WFIRST_CPU=0 ./scripts/benchmark_all.sh

# Heavy WFIRST parity matrix.
PYTHON_BIN=.venv-parity/bin/python WFIRST_PARITY_ONLY=1 ./scripts/benchmark_wfirst_phaseb_cpu.sh --parity-only
```

Julia benchmark dependencies live in `bench/Project.toml`, not the core package
environment. The benchmark scripts automatically develop the local checkout into
that environment before running Julia benchmark code.

Optional GPU benchmark backends should be added to the benchmark environment,
not to the core package environment:

```bash
julia --project=bench -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add("CUDA")'
julia --project=bench -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add("AMDGPU")'
```

For a targeted WFIRST check:

```bash
PYTHON_BIN=.venv-parity/bin/python \
WFIRST_PARITY_ONLY=1 \
WFIRST_CASES=full_none,full_hlc \
./scripts/benchmark_wfirst_phaseb_cpu.sh --parity-only
```

### Julia CPU vs GPU Benchmarks
Run the dedicated Julia CPU/GPU comparison lane with:

```bash
./scripts/benchmark_cpu_gpu.sh
```

This generates:
- `bench/reports/julia_cpu_gpu_summary.md`
- `bench/reports/julia_cpu_gpu_steady_state.csv`
- `bench/reports/julia_cpu_gpu_batch_throughput.csv`
- `bench/reports/julia_cpu_gpu_precision_split.csv`
- `bench/reports/julia_cpu_gpu_core_propagation_tail.csv`
- `bench/reports/julia_cpu_gpu_supported_kernels.csv`

The script uses the Julia steady-state CPU lane plus any available GPU lanes.
It does not depend on the Python parity environment. If a GPU backend or
supported device is unavailable, it records that backend as skipped instead of
failing the whole run. GPU benchmark lanes require `CUDA.jl` or `AMDGPU.jl` to
be available in the `bench/` environment on the machine running the benchmark.

The summary includes both:
- the standard steady-state workload
- prepared batch throughput for repeated prepared runs
- FP64/FP32 precision comparisons for supported backends
- a synthetic core propagation-tail workload derived from the repeated
  lens/propagate sequence identified during hotspot analysis

Representative validated batch result on CUDA hardware:
- 4-wavelength prepared sweep, `512 x 512` grid
- CPU FP64: `103.42 ms`
- CUDA FP64: `7.82 ms`
- CPU/CUDA FP64: `13.22x`
- CPU FP32: `70.11 ms`
- CUDA FP32: `537.67 us`
- CPU/CUDA FP32: `130.40x`

That is the clearest current example of why prepared batched execution and
explicit `Float32` support are first-class performance features rather than
benchmark-only details.

For profiling the shared core workload directly, run:

```bash
./scripts/profile_core_cpu_gpu.sh
```

This writes backend-specific text profiles under `bench/reports/` and keeps GPU
optimization work tied to shared propagation code rather than model-specific
wrappers.

For Python-baseline parity and Python-vs-Julia benchmark runs, install the
official Python PROPER 3.3.4 source tree first. For WFIRST/Roman Phase B
Python-vs-Julia comparisons, also install the upstream `proper-models` source
tree:

```bash
./scripts/setup_python_baseline.sh
./scripts/setup_wfirst_models_baseline.sh
./scripts/setup_parity_venv.sh
PYTHON_BIN=.venv-parity/bin/python ./scripts/benchmark_all.sh
```

Set `PYPROPER_ROOT=/path/to/proper_v3.3.4_python` if the baseline lives outside
the default sibling checkout path. Set `WFIRST_MODELS_ROOT=/path/to/proper-models`
if the WFIRST/Roman model checkout lives outside the default sibling path. The
WFIRST Phase B Python source does not require Git LFS assets; for a full
`proper-models` checkout with LFS assets, run
`WFIRST_MODELS_LFS=1 ./scripts/setup_wfirst_models_baseline.sh`.

### GPU Usage Contract
- The intended GPU performance surface is:
  - mutating `!` APIs
  - explicit `RunContext`
  - prepared execution for repeated runs
- Allocating wrappers remain available for convenience, but they are not the
  steady-state performance contract.
- `BenchmarkTools` results in the benchmark lanes are warmed steady-state
  numbers, not cold-start / TTFx numbers.

Examples:

```julia
ctx = RunContext(typeof(wf.field))
prop_resamplemap!(out, wf, dmap, opts, ctx)
prop_magnify!(out, image_in, mag, ctx; QUICK=true)
```

Current GPU support is strongest for:
- propagation core (`prop_qphase`, `prop_ptp`, `prop_wts`, `prop_stw`, `prop_end!`)
- same-backend map application (`prop_add_phase`, `prop_multiply`, `prop_divide`)
- promoted map/error-map paths (`prop_readmap`, `prop_errormap`, `prop_psd_errormap`)

Unsupported GPU combinations fail explicitly instead of silently falling back to
host materialization. See [docs/backend_traits.md](docs/backend_traits.md) for
the current support matrix.

## Requirements
- Julia 1.10 or newer
- `FITSIO.jl` for FITS input/output

Package dependencies are declared in `Project.toml` and resolved through Julia's
package manager.

## Installation

### Use As A Project Dependency
From Julia:

```julia
import Pkg
Pkg.add(path="/path/to/proper.jl")
```

Or for development:

```julia
import Pkg
Pkg.develop(path="/path/to/proper.jl")
```

### Use Without Installing Globally
Activate the project and run Julia with `--project=.`:

```bash
julia --project=.
```

Then in Julia:

```julia
using Proper
```

## Plotting Output

Plotting is intentionally not part of the core `Proper.jl` runtime dependency
set. Example scripts use `Plots.jl` through `examples/Project.toml`.

For local examples from a checkout:

```bash
julia --project=examples -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=examples examples/simple_prescription.jl
```

For your own application code, choose the plotting stack you prefer. With
`Plots.jl`, a typical PSF display is:

```julia
using Plots
heatmap(log10.(abs.(psf) .+ eps()); aspect_ratio=:equal, title="PSF")
```

## Writing FITS Output

```julia
prop_fits_write("example.fits", psf)
```

Map-oriented output can be written with PROPER-compatible metadata using
`prop_writemap`.

## Documentation
- [Documentation index](docs/README.md)
- [API contract](docs/api_contract.md)
- [Migration guide](docs/MIGRATION_GUIDE.md)
- [Prescription authoring guide](docs/PRESCRIPTION_AUTHORING_GUIDE.md)
- [Contributor porting checklist](docs/PORTING_CHECKLIST.md)
- [Prepared execution guide](docs/PREPARED_EXECUTION_GUIDE.md)
- [Runnable API examples](docs/API_EXAMPLES.md)
- [Numerics contract](docs/numerics_contract.md)
- [Compatibility decisions](docs/compat_decisions.md)
