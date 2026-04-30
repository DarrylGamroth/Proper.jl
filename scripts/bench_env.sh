#!/usr/bin/env bash

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROPER_REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
PROPER_BENCH_PROJECT="${PROPER_REPO_ROOT}/bench"

ensure_bench_env() {
  if [[ "${PROPER_BENCH_ENV_READY:-0}" == "1" ]]; then
    return 0
  fi

  if [[ "${PROPER_BENCH_SKIP_INSTANTIATE:-0}" != "1" ]]; then
    PROPER_REPO_ROOT="${PROPER_REPO_ROOT}" julia --project="${PROPER_BENCH_PROJECT}" -e '
      using Pkg
      Pkg.develop(Pkg.PackageSpec(path=ENV["PROPER_REPO_ROOT"]))
      Pkg.instantiate()
    '
  fi

  export PROPER_BENCH_ENV_READY=1
}

bench_julia() {
  ensure_bench_env
  julia --project="${PROPER_BENCH_PROJECT}" "$@"
}
