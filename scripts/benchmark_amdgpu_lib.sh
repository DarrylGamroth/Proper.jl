#!/usr/bin/env bash

run_amdgpu_benchmarks() {
  local probe_output=""

  if probe_output=$(julia --project=. bench/julia/amdgpu/probe.jl 2>&1); then
    run_step "Julia AMDGPU steady-state workload" julia --project=. bench/julia/amdgpu/steady_state.jl
    run_step "Julia AMDGPU core propagation tail" julia --project=. bench/julia/amdgpu/core_propagation_tail.jl
    run_step "Julia AMDGPU supported kernels" julia --project=. bench/julia/amdgpu/supported_kernels.jl
  else
    echo "[bench] AMDGPU benchmarks skipped"
    AMDGPU_SKIP_REASON="${probe_output}" julia --project=. bench/julia/amdgpu/write_skipped_reports.jl >/dev/null
  fi
}
