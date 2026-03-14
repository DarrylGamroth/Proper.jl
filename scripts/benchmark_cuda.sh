#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/benchmark_cuda_lib.sh"

run_step() {
  local label="$1"
  shift

  echo "[bench] ${label}"
  if [[ "${BENCH_VERBOSE:-0}" == "1" ]]; then
    "$@"
  else
    "$@" >/dev/null
  fi
}

run_cuda_benchmarks
run_step "Julia cold-start / TTFx" julia --project=. bench/julia/cold_start/run.jl
julia --project=. bench/reports/summarize.jl
