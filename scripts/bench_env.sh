#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROPER_REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
PROPER_BENCH_PROJECT="${PROPER_BENCH_PROJECT:-${PROPER_REPO_ROOT}/bench}"
PROPER_LATENCY_PROJECT="${PROPER_LATENCY_PROJECT:-${PROPER_REPO_ROOT}/bench/latency}"

ensure_bench_env() {
  if [[ "${PROPER_BENCH_ENV_READY:-0}" == "1" ]]; then
    return 0
  fi

  if [[ "${PROPER_BENCH_SKIP_INSTANTIATE:-0}" != "1" ]]; then
    PROPER_REPO_ROOT="${PROPER_REPO_ROOT}" julia --project="${PROPER_BENCH_PROJECT}" -e '
      using Pkg
      using TOML
      project = TOML.parsefile(Base.active_project())
      sources = get(project, "sources", Dict{String,Any}())
      if VERSION < v"1.11" || !haskey(sources, "Proper")
        Pkg.develop(Pkg.PackageSpec(path=ENV["PROPER_REPO_ROOT"]))
      end
      Pkg.instantiate()
    '
  fi

  export PROPER_BENCH_ENV_READY=1
}

bench_julia() {
  ensure_bench_env
  julia --project="${PROPER_BENCH_PROJECT}" "$@"
}

ensure_latency_env() {
  ensure_bench_env
  if [[ "${PROPER_LATENCY_ENV_READY:-0}" == "1" ]]; then
    return 0
  fi

  if [[ "${PROPER_BENCH_SKIP_INSTANTIATE:-0}" != "1" ]]; then
    julia --startup-file=no --project="${PROPER_LATENCY_PROJECT}" -e '
      using Pkg
      Pkg.instantiate()
    '
  fi

  export PROPER_LATENCY_ENV_READY=1
}

latency_julia() {
  ensure_latency_env
  JULIA_LOAD_PATH="${PROPER_BENCH_PROJECT}:${PROPER_LATENCY_PROJECT}:@stdlib" \
    julia --startup-file=no --project="${PROPER_BENCH_PROJECT}" "$@"
}
