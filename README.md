# Proper.jl

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
- Phase 8 parity closure is complete against the patched Python 3.3.4 baseline:
  see [docs/PHASE8_CLOSURE.md](docs/PHASE8_CLOSURE.md)
- Phase 9 MATLAB/manual semantic reconciliation is complete for the identified
  hotspot semantics: see [docs/PHASE9_RECONCILIATION.md](docs/PHASE9_RECONCILIATION.md)
- The WFIRST Phase B reference port is used as a broad correctness and
  benchmarking workload, not as a separate optimized execution path
- Current benchmark policy separates steady-state runtime from Julia cold-start /
  TTFx (`D-0029`)

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

### Use `prepare_prescription_batch(...)` For Repeated Parallel Work
Use `PreparedBatch` when you want reusable per-slot contexts for repeated or
parallel execution.

```julia
batch = prepare_prescription_batch(simple_prescription, 0.55, 128; pool_size=2)
stack, samplings = prop_run_multi(batch; PASSVALUE=[nothing, nothing])
```

### Use `prepare_model(...)` When Assets Or Naming Matter
Use `PreparedModel` when you want a named execution object and optional cached
assets layered on top of prepared contexts.

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

## Documentation
- [API contract](docs/api_contract.md)
- [Migration guide](docs/MIGRATION_GUIDE.md)
- [Prescription authoring guide](docs/PRESCRIPTION_AUTHORING_GUIDE.md)
- [Contributor porting checklist](docs/PORTING_CHECKLIST.md)
- [Prepared execution guide](docs/PREPARED_EXECUTION_GUIDE.md)
- [Runnable API examples](docs/API_EXAMPLES.md)
- [Numerics contract](docs/numerics_contract.md)
- [Compatibility decisions](docs/compat_decisions.md)
