#!/usr/bin/env bash

run_amdgpu_benchmarks() {
  local probe_output=""

  if probe_output=$(bench_julia bench/julia/amdgpu/probe.jl 2>&1); then
    run_step "Julia AMDGPU steady-state workload" bench_julia bench/julia/amdgpu/steady_state.jl
    run_step "Julia AMDGPU steady-state FP32" bench_julia bench/julia/amdgpu/steady_state_fp32.jl
    run_step "Julia AMDGPU batch throughput FP64" bench_julia bench/julia/amdgpu/batch_throughput_fp64.jl
    run_step "Julia AMDGPU batch throughput FP32" bench_julia bench/julia/amdgpu/batch_throughput_fp32.jl
    run_step "Julia AMDGPU core propagation tail" bench_julia bench/julia/amdgpu/core_propagation_tail.jl
    if run_step "Julia AMDGPU supported kernels" bench_julia bench/julia/amdgpu/supported_kernels.jl; then
      :
    else
      echo "[bench] AMDGPU supported-kernel lane skipped after compiler failure"
      AMDGPU_SKIP_REASON="AMDGPU helper-kernel benchmark failed under the current compiler/toolchain" \
        bench_julia bench/julia/amdgpu/write_skipped_supported_kernels.jl >/dev/null
    fi
  else
    echo "[bench] AMDGPU benchmarks skipped"
    AMDGPU_SKIP_REASON="${probe_output}" bench_julia bench/julia/amdgpu/write_skipped_reports.jl >/dev/null
  fi
}
