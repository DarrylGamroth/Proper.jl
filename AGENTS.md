# AGENTS.md

## Scope
This file defines local agent/contributor guidance for `proper.jl` (PROPER Python/MATLAB -> Julia port).

## Source Of Truth
- Port plan: `PORTING_PLAN.md`
- Compatibility decisions: `docs/compat_decisions.md`
- Python baseline code: `../proper_v3.3.4_python`
- MATLAB semantic reference: `../proper_v3.3.1_matlab` (no runtime available here)

## Python Parity Environment
- Use local venv: `.venv-parity`.
- Baseline/parity generators should run with `.venv-parity/bin/python`.
- Required Python deps for parity workflows: `numpy`, `scipy`, `astropy`, `matplotlib`.

## Non-Negotiable Decisions (Accepted)
- `D-0001`: Python 3.3.4 is the executable parity baseline.
- `D-0002`: FITS handling uses `FITSIO.jl`.
- `D-0003`: Default plotting for examples uses `Plots.jl`.
- `D-0029`: Benchmark reporting must separate steady-state runtime from Julia cold-start/TTFx.

## Porting Policy
- Keep one-to-one filename mapping with Python for traceability.
- Do not do line-by-line transliteration; use idiomatic Julia internals.
- Preserve familiar public `prop_*` API names for PROPER users.
- Prefer multiple dispatch, concrete types, and backend traits.
- Default behavior should remain parity-first unless an explicit decision says otherwise.
- Hot-path code must be type-stable and minimize dynamic dispatch.

## Compatibility Baseline
- Behavior targets the patched Python 3.3.4 executable baseline used by parity harnesses.
- MATLAB/manual remain semantic references when evaluating suspected translation defects.
- Do not add runtime compatibility mode flags; document and test behavior changes directly.

## Array/Backend Requirements
- Core APIs should accept `AbstractArray`/`AbstractMatrix`.
- Preserve input backend (`Array`, `CuArray`, etc.) where feasible.
- Route backend-specific behavior via traits and dispatch, not mutable global flags.
- Use `KernelAbstractions.jl`/`AcceleratedKernels.jl` where appropriate for portable kernels.

## Contribution Rules For Agents
- Before coding:
  - Check `PORTING_PLAN.md` Phase 0/phase status.
  - Check `docs/compat_decisions.md` for accepted/proposed decisions.
- When behavior differs from Python:
  - Add/update a decision entry with context, choice, and consequences.
  - Add tests that cover both compatibility and corrected behavior (if applicable).
- When porting modules/examples:
  - Keep naming/folder parity.
  - Add or update parity tests and provenance metadata.
  - Validate type stability for hot kernels and avoid runtime symbol/dict branching inside loops.

## Testing Expectations
- Regression parity compares Julia outputs against Python baselines.
- Randomized paths must be seeded during tests/parity generation.
- GPU tests must avoid scalar indexing (`CUDA.allowscalar(false)` when applicable).
