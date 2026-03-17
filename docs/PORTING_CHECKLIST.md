# Porting Checklist

## Purpose
This checklist is for contributors porting new PROPER modules, examples, or
model code into `Proper.jl`.

It turns the project policy in `PORTING_PLAN.md` and
`docs/compat_decisions.md` into a concrete workflow.

## Before You Start
- read the relevant phase scope in `PORTING_PLAN.md`
- check `docs/compat_decisions.md` for accepted behavior choices
- identify the upstream source file in:
  - `../proper_v3.3.4_python`
- identify whether MATLAB source in:
  - `../proper_v3.3.1_matlab`
  is likely to matter for semantic review

## File And Naming Rules
- keep one-to-one filename mapping with the upstream Python source
- preserve public `prop_*` naming for user-facing routines
- do not do line-by-line transliteration if Julia has a clearer typed design
- keep traceability obvious in file naming and test naming

## Implementation Checklist

### 1. Establish The Boundary
- define the public entry point with compatibility-facing `prop_*` naming
- accept uppercase compatibility keywords where required
- normalize options once at the API boundary
- keep hot-path internals typed and dispatch-driven

### 2. Keep Internals Julia-Native
- replace upstream module-global mutable state with explicit arguments or context
- prefer concrete structs, `NamedTuple`s, and multiple dispatch
- avoid `Any`-typed state in hot paths
- avoid runtime string/symbol branching inside loops

### 3. Preserve Backend Flexibility
- accept `AbstractArray` / `AbstractMatrix` where practical
- do not force `Array` unless host materialization is actually required
- route backend-specific behavior through traits and dispatch
- avoid GPU-hostile scalar indexing in reusable kernels

### 4. Make Semantics Explicit
- check centering, indexing, and unit conventions
- if Python and MATLAB disagree, determine whether this is:
  - executable-baseline behavior to preserve, or
  - a translation defect to correct
- document any accepted divergence in `docs/compat_decisions.md`

### 5. Prefer Reusable Workspaces
- use existing `RunContext` / workspace facilities for hot paths
- avoid per-call temporary allocations in steady-state kernels
- use mutating `...!` internals when they improve ownership and performance

## Testing Checklist

### Required
- add or update unit tests for the ported routine
- add parity coverage if the routine participates in parity workflows
- seed randomized paths in tests
- keep return shape and keyword behavior covered

### When Behavior Differs From Python
- add a decision entry in `docs/compat_decisions.md`
- add tests that lock in the accepted behavior
- reference MATLAB/manual semantics when that is the reason for the difference

### When Performance Matters
- check type stability for hot-path code
- add allocation/inference guards where appropriate
- benchmark steady-state behavior separately from Julia cold-start / TTFx

## Parity Workflow
- use `.venv-parity/bin/python` for executable Python baseline runs
- keep Python 3.3.4 as the primary executable parity reference
- use MATLAB as the semantic reference for suspected translation defects,
  especially column-major-sensitive routines

## Prescription And Example Ports
- port the optical sequence first using familiar public calls
- get plain `prop_run(...)` working before introducing prepared execution
- use prepared execution only after parity is established
- keep examples readable; they are user-facing documentation as well as tests

## FITS And External Data
- use the established FITS helpers
- keep array-order policy explicit
- localize compatibility reshaping to the model layer when possible
- do not silently bake model-specific FITS assumptions into shared core code

## Review Checklist
- file mapping preserved
- public API remains familiar
- behavior choice documented if it differs from Python
- tests added or updated
- parity evidence updated when relevant
- docs updated if the new port changes user-visible behavior

## Related Docs
- `PORTING_PLAN.md`
- `docs/compat_decisions.md`
- `docs/api_contract.md`
- `docs/MIGRATION_GUIDE.md`
- `docs/PRESCRIPTION_AUTHORING_GUIDE.md`
