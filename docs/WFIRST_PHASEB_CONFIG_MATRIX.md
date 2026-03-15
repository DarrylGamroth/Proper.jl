# WFIRST Phase B Config Matrix

Representative configuration-matrix coverage for the Julia reference port and the Python parity harness.

| Area | Axis | Permutations | Evidence | Status | Notes |
| --- | --- | --- | --- | --- | --- |
| Coronagraph family | `cor_type` | `hlc`, `spc-spec_long`, `spc-wide`, `none` | `scripts/verify_wfirst_phaseb_matrix.sh` | Covered | `hlc_erkin`, `spc-ifs_short`, and `spc-ifs_long` still need parity cases. |
| Model size | prescription path | `compact`, `full` | `scripts/verify_wfirst_phaseb_matrix.sh` | Covered | Matrix includes both compact and full runs for HLC and SPC. |
| Source offset | source tilt | zero, nonzero lambda-D | `compact_hlc_source_offset` | Gap | Current compact HLC parity case exposes a real Python-vs-Julia mismatch and remains under investigation. |
| Source offset conversion | mas to `lambda/D` | zero, nonzero mas | `test/test_wfirst_phaseb_reference.jl` | Covered | Helper-level coverage; not yet a separate Python-vs-Julia parity row. |
| HLC field stop | `use_field_stop` | `1`, `0` | `full_hlc`, `full_hlc_no_field_stop` | Covered | Full HLC path only. |
| SPC pupil mask | `use_pupil_mask` | `1`, `0` | `full_spc_spec_long`, `full_spc_spec_long_no_pupil_mask` | Covered | Spec-long branch only. |
| Coronagraph bypass | `use_fpm` / `none` | normal coronagraph, bypassed optics | `full_none` | Covered | `full_none` exercises the pupil-only path. |
| Error maps | `use_errors` | `0`, `1` | reference tests + CPU parity harness | Gap | Current parity surface validates `use_errors=0` only. |
| DMs | `use_dm1`, `use_dm2`, explicit maps | off, on with supplied maps | unit coverage in `test/test_wfirst_phaseb_reference.jl` | Gap | Helper/argument handling is covered; full Python-vs-Julia parity cases are still missing. |
| SPC family variants | SPC branch | `spc-ifs_short`, `spc-ifs_long` | none yet | Gap | Data/config plumbing exists, parity cases not yet added. |
| HLC variant | `hlc_erkin` | default, alt HLC branch | none yet | Gap | Config selection is implemented but not yet exercised in parity harness. |

## Standard Verification

Run the representative parity matrix without timing loops:

```bash
JULIA_NUM_THREADS=4 ./scripts/verify_wfirst_phaseb_matrix.sh
```

Run the full WFIRST CPU benchmark with timings:

```bash
JULIA_NUM_THREADS=4 ./scripts/benchmark_wfirst_phaseb_cpu.sh
```

Select a subset or parity-only mode directly:

```bash
JULIA_NUM_THREADS=4 ./scripts/benchmark_wfirst_phaseb_cpu.sh \
  --cases compact_hlc,compact_hlc_source_offset,full_none \
  --parity-only
```
