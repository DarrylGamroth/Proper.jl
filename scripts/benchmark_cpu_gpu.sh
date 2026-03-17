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

run_step "Julia CPU steady-state workload" julia --project=. bench/julia/steady_state/run.jl
run_step "Julia CPU supported kernels" julia --project=. bench/julia/steady_state/supported_kernels.jl
run_cuda_benchmarks
julia --project=. bench/reports/summarize_cpu_gpu.jl
