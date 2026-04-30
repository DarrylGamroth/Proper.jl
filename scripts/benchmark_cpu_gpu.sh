#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/bench_env.sh"
source "$(dirname "$0")/benchmark_cuda_lib.sh"
source "$(dirname "$0")/benchmark_amdgpu_lib.sh"

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

run_step "Julia CPU steady-state workload" bench_julia bench/julia/steady_state/run.jl
run_step "Julia CPU steady-state FP32" bench_julia bench/julia/steady_state/run_fp32.jl
run_step "Julia CPU batch throughput FP64" bench_julia bench/julia/steady_state/batch_throughput_fp64.jl
run_step "Julia CPU batch throughput FP32" bench_julia bench/julia/steady_state/batch_throughput_fp32.jl
run_step "Julia CPU core propagation tail" bench_julia bench/julia/steady_state/core_propagation_tail.jl
run_step "Julia CPU supported kernels" bench_julia bench/julia/steady_state/supported_kernels.jl
run_cuda_benchmarks
run_amdgpu_benchmarks
bench_julia bench/reports/summarize_cpu_gpu.jl
