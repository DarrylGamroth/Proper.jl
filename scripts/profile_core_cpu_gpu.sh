#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/benchmark_cuda_lib.sh"
source "$(dirname "$0")/benchmark_amdgpu_lib.sh"

profile_step() {
  local label="$1"
  shift

  echo "[profile] ${label}"
  if [[ "${BENCH_VERBOSE:-0}" == "1" ]]; then
    "$@"
  else
    "$@" >/dev/null
  fi
}

profile_step "Julia CPU core propagation tail" julia --project=. bench/julia/core/profile_propagation_tail.jl cpu

if julia --project=. bench/julia/cuda/probe.jl >/dev/null 2>&1; then
  profile_step "Julia CUDA core propagation tail" julia --project=. bench/julia/core/profile_propagation_tail.jl cuda
else
  echo "[profile] CUDA core propagation tail skipped"
fi

if julia --project=. bench/julia/amdgpu/probe.jl >/dev/null 2>&1; then
  profile_step "Julia AMDGPU core propagation tail" julia --project=. bench/julia/core/profile_propagation_tail.jl amdgpu
else
  echo "[profile] AMDGPU core propagation tail skipped"
fi
