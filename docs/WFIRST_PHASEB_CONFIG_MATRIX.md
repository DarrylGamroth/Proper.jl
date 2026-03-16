# WFIRST Phase B Config Matrix

Representative configuration-matrix coverage for the Julia reference port and the Python parity harness.

| Area | Axis | Permutations | Evidence | Status | Notes |
| --- | --- | --- | --- | --- | --- |
| Coronagraph family | `cor_type` | `hlc`, `spc-spec_short`, `spc-spec_long`, `spc-ifs_short`, `spc-ifs_long`, `spc-wide`, `none` | `scripts/verify_wfirst_phaseb_matrix.sh` | Covered | `hlc_erkin` still needs dedicated parity cases and compatible public data. |
| Model size | prescription path | `compact`, `full` | `scripts/verify_wfirst_phaseb_matrix.sh` | Covered | Matrix includes both compact and full runs for HLC and SPC. |
| Source offset | source tilt | zero, nonzero lambda-D | `compact_hlc_source_offset`, `full_hlc_source_offset` | Covered | HLC source-offset parity now matches at machine precision on both compact and full rows. |
| Source offset conversion | mas to `lambda/D` | zero, nonzero mas | `test/test_wfirst_phaseb_reference.jl` | Covered | Helper-level coverage; not yet a separate Python-vs-Julia parity row. |
| HLC field stop | `use_field_stop` | `1`, `0` | `full_hlc`, `full_hlc_no_field_stop` | Covered | Full HLC path only. |
| SPC pupil mask | `use_pupil_mask` | `1`, `0` | `full_spc_spec_long`, `full_spc_spec_long_no_pupil_mask` | Covered | Spec-long branch only. |
| Coronagraph bypass | `use_fpm` / `none` | normal coronagraph, bypassed optics | `full_none` | Covered | `full_none` exercises the pupil-only path. |
| Error maps | `use_errors` | `0`, `1` | reference tests + CPU parity harness | Gap | Current parity surface validates `use_errors=0` only. |
| DMs | `use_dm1`, `use_dm2`, explicit maps | off, on with supplied maps | `compact_hlc_dm_pair`, `full_hlc_dm_pair` | Covered | Explicit HLC DM-map parity rows now match with relative L2 around `4.13e-3`; this is treated as comparable fidelity rather than exact internal equivalence. |
| SPC family variants | SPC branch | `spc-spec_short`, `spc-spec_long`, `spc-ifs_short`, `spc-ifs_long`, `spc-wide` | `scripts/verify_wfirst_phaseb_matrix.sh` | Covered | Current public-data matrix validates both spec and IFS aliases. |
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
