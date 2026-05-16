# Performance Card Stack

This directory holds card-specific benchmark helpers, prompt notes, and other
temporary/per-card assets for the performance card stack. Card-only helpers
must stay here or in a card-specific subdirectory such as `bench/cards/card00/`.
Do not put card helpers in shared source, package, or general benchmark code
directories unless they have become reusable infrastructure.

## Fresh Codex Context

Use a fresh Codex/ChatGPT context for each performance card. Shared state should
come from git history and files in the repository, not from prior chat memory.

For each card prompt:

1. Start a new Codex/ChatGPT task for that card.
2. Tell the agent to re-read `AGENTS.md`, `PORTING_PLAN.md`, and
   `docs/compat_decisions.md`.
3. Start from clean updated `main`.
4. Create only that card's branch.
5. Run baseline benchmarks on `main` before implementation.
6. Keep the card independently reviewable, revertible, and mergeable.
7. Save benchmark reports under ignored `bench/reports/card_<id>_*.json`.

If the temporary Card 00 harness is still uncommitted when starting a new card,
stash it before updating `main`, create the new card branch from clean `main`,
then re-apply the stash on that branch:

```bash
git stash push -u -m "temporary card00 benchmark harness"
git switch main
git pull --ff-only
git status --short
git switch -c perf/card-XX-short-name
git stash apply stash^{/"temporary card00 benchmark harness"}
```

## Branch Template

```bash
git switch main
git pull --ff-only
git status --short
git switch -c perf/card-XX-short-name
```

Only merge a card when parity/tests pass and the selected gate is met: at least
10% median steady-state runtime improvement or 25% allocation reduction.

## Card 00 Layout

Card 00 benchmark helpers live in `bench/cards/card00/`.

Benchmark entry points currently live with the existing benchmark lanes:

- `bench/julia/steady_state/card00_dm_projection.jl`
- `bench/julia/steady_state/card00_zernike_synthesis.jl`
- `bench/julia/steady_state/card00_zernike_fit.jl`
- `bench/julia/wfirst_phaseb/card00_prepared_models.jl`
