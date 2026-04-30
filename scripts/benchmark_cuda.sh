#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROPER_REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
export PROPER_BENCH_PROJECT="${PROPER_BENCH_PROJECT:-${PROPER_REPO_ROOT}/bench/cuda}"

source "${SCRIPT_DIR}/bench_env.sh"
source "${SCRIPT_DIR}/benchmark_cuda_lib.sh"

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
