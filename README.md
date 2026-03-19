# Proper.jl

[![CI](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/DarrylGamroth/Proper.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/DarrylGamroth/Proper.jl/graph/badge.svg)](https://codecov.io/gh/DarrylGamroth/Proper.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia port of John Krist's PROPER optical propagation library.

PROPER is a diffraction propagation toolkit for exploring wave-optics behavior
in optical systems. It is not a ray-tracing system and is not intended to
replace detailed optical design tools. `Proper.jl` keeps the familiar public
`prop_*` API while using Julia-native internals, explicit runtime state, and a
prepared execution layer for repeated runs.

## Project Status
- public `prop_*` routine naming preserved for migration familiarity
- Python PROPER 3.3.4 is the executable parity baseline
- MATLAB PROPER 3.3.1 is used as a semantic reference, especially for
  column-major correctness questions
- broad CPU parity coverage is in place, including WFIRST Phase B reference
  workloads
- prepared execution (`PreparedPrescription`, `PreparedBatch`,
  `PreparedModel`) is available for repeated and parallel runs

## Benchmark And Parity Summary
- Parity closure is complete against the patched Python 3.3.4 baseline:
  see [docs/PHASE8_CLOSURE.md](docs/PHASE8_CLOSURE.md)
- MATLAB/manual semantic reconciliation is complete for the identified hotspot
  semantics: see [docs/PHASE9_RECONCILIATION.md](docs/PHASE9_RECONCILIATION.md)
- The WFIRST Phase B reference port is used as a broad correctness and
  benchmarking workload, not as a separate optimized execution path
- Current benchmark policy separates steady-state runtime from Julia cold-start /
  TTFx (`D-0029`)

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

The script uses the Julia steady-state CPU lane plus any available GPU lanes
(`CUDA.jl` and/or `AMDGPU.jl`). It does not depend on the Python parity
environment. If a GPU backend or supported device is unavailable, it records
that backend as skipped instead of failing the whole run.

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
- `Plots.jl` for example plotting

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

## Interactive Use

The usual migration path is:

1. write or port a prescription using the familiar `prop_*` calls
2. run it with `prop_run(...)`
3. compare against the Python baseline if parity matters

You can also run multiple cases:

```julia
stack, samplings = prop_run_multi(simple_prescription, 0.55, 128; PASSVALUE=[nothing, nothing])
```

`prop_run_multi` preserves input order and returns a 3-D output stack plus a
vector of samplings.

## Plotting Output

Example plotting uses `Plots.jl` by default:

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

## Choose The Right Execution Surface

### Use `prop_run(...)` First
Use plain `prop_run` when:
- porting or validating a prescription
- comparing behavior against upstream Python
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
- the WFIRST Phase B reference port is included as a broad correctness and
  benchmarking workload
- the WFIRST reference model lives under `reference_models/wfirst_phaseb_proper`

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
