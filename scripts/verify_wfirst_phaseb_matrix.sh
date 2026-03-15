#!/usr/bin/env bash
set -euo pipefail

cases="${WFIRST_MATRIX_CASES:-compact_hlc,full_hlc,compact_spc_spec_long,full_spc_spec_long,compact_spc_wide,full_spc_wide,full_hlc_no_field_stop,full_spc_spec_long_no_pupil_mask,full_none}"

exec ./scripts/benchmark_wfirst_phaseb_cpu.sh \
  --cases "${cases}" \
  --samples 0 \
  --parity-only
