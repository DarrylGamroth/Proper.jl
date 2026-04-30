#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/bench_env.sh"
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
run_step "Julia cold-start / TTFx" bench_julia bench/julia/cold_start/run.jl
bench_julia bench/reports/summarize.jl
