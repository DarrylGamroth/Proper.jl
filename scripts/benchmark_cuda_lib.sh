#!/usr/bin/env bash

run_cuda_benchmarks() {
  local probe_output=""

  if probe_output=$(julia --project=. bench/julia/cuda/probe.jl 2>&1); then
    run_step "Julia CUDA steady-state workload" julia --project=. bench/julia/cuda/steady_state_fp64.jl
    run_step "Julia CUDA core propagation tail" julia --project=. bench/julia/cuda/core_propagation_tail.jl
    run_step "Julia CUDA steady-state FP32" julia --project=. bench/julia/cuda/steady_state_fp32.jl
    run_step "Julia CUDA batch throughput FP64" julia --project=. bench/julia/cuda/batch_throughput_fp64.jl
    run_step "Julia CUDA batch throughput FP32" julia --project=. bench/julia/cuda/batch_throughput_fp32.jl
    run_step "Julia CUDA supported kernels" julia --project=. bench/julia/cuda/supported_kernels.jl
    run_step "Julia CUDA precision split" julia --project=. bench/julia/cuda/precision_split.jl

    local kernels=(
      prop_qphase
      prop_ptp
      prop_wts
      prop_stw
      prop_circular_aperture
      prop_end_mutating
    )
    local precision=""
    local kernel=""
    for precision in fp64 fp32; do
      for kernel in "${kernels[@]}"; do
        run_step "Julia CUDA isolated ${precision} ${kernel}" julia --project=. bench/julia/cuda/isolated_wavefront_kernel.jl "${kernel}" "${precision}"
      done
    done

    run_step "Julia CUDA isolated wavefront summary" julia --project=. bench/julia/cuda/aggregate_isolated_wavefront_kernels.jl
  else
    echo "[bench] CUDA benchmarks skipped"
    CUDA_SKIP_REASON="${probe_output}" julia --project=. bench/julia/cuda/write_skipped_reports.jl >/dev/null
  fi
}
