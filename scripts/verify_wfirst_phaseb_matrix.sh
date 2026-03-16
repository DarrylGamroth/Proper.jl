#!/usr/bin/env bash
set -euo pipefail

cases="${WFIRST_MATRIX_CASES:-compact_hlc,full_hlc,full_hlc_errors,compact_hlc_source_offset,full_hlc_source_offset,compact_hlc_dm_pair,full_hlc_dm_pair,compact_spc_spec_short,full_spc_spec_short,compact_spc_ifs_short,full_spc_ifs_short,compact_spc_spec_long,full_spc_spec_long,full_spc_spec_long_errors,compact_spc_ifs_long,full_spc_ifs_long,compact_spc_wide,full_spc_wide,full_hlc_no_field_stop,full_spc_spec_long_no_pupil_mask,full_none}"

exec ./scripts/benchmark_wfirst_phaseb_cpu.sh \
  --cases "${cases}" \
  --samples 0 \
  --parity-only
