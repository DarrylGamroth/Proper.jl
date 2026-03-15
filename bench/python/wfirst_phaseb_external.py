#!/usr/bin/env python3
import argparse
import json
import os
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def write_skipped_report(case_name: str, reason: str) -> None:
    outdir = project_root() / "bench" / "reports"
    outdir.mkdir(parents=True, exist_ok=True)
    report = {
        "meta": {
            "run_tag": f"wfirst_phaseb_python_{case_name}",
            "backend": "cpu",
            "baseline": "wfirst_phaseb_python",
            "case": case_name,
            "available": False,
        },
        "reason": reason,
    }
    outpath = outdir / f"python_wfirst_phaseb_{case_name}.json"
    outpath.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


def resolve_roots() -> tuple[Path, Path, Path | None]:
    root = project_root()
    proper_root = (root / ".." / "proper_v3.3.4_python").resolve()
    models_python_root = (root / ".." / "proper-models" / "wfirst_cgi" / "models_phaseb" / "python").resolve()
    env_data_root = os.environ.get("WFIRST_PHASEB_DATA_ROOT")
    candidate_data_roots = []
    if env_data_root:
        candidate_data_roots.append(Path(env_data_root).expanduser().resolve())
    candidate_data_roots.extend([
        (root / ".." / "proper-models" / "wfirst_cgi" / "data_phaseb").resolve(),
        (root / ".." / "proper-models" / "wfirst_cgi" / "models_phaseb" / "data_phaseb").resolve(),
    ])
    data_root = next((path for path in candidate_data_roots if path.is_dir()), None)

    if not proper_root.is_dir():
        raise RuntimeError(f"Missing Python PROPER baseline source tree: {proper_root}")
    if not models_python_root.is_dir():
        raise RuntimeError(f"Missing WFIRST Phase B Python source tree: {models_python_root}")
    return proper_root, models_python_root, data_root


def load_python_models():
    proper_root, models_python_root, data_root = resolve_roots()
    sys.path.insert(0, str(proper_root))
    sys.path.insert(0, str(models_python_root))

    try:
        import proper  # noqa: WPS433
        import wfirst_phaseb_proper  # noqa: WPS433
        from wfirst_phaseb_proper.wfirst_phaseb import wfirst_phaseb  # noqa: WPS433
        from wfirst_phaseb_proper.wfirst_phaseb_compact import wfirst_phaseb_compact  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - external harness guardrail
        raise RuntimeError(
            "Unable to import Python PROPER or WFIRST Phase B model. "
            f"Interpreter: {sys.executable}"
        ) from exc

    proper.print_it = False
    proper.verbose = False
    proper.print_total_intensity = False
    proper.do_table = False

    if data_root is not None:
        wfirst_phaseb_proper.data_dir = str(data_root)

    return proper, wfirst_phaseb, wfirst_phaseb_compact, proper_root, models_python_root, data_root


def case_definitions(wfirst_phaseb, wfirst_phaseb_compact):
    hlc_lam0_um = 0.575
    hlc_band = 0.1
    hlc_lams_um = np.linspace(hlc_lam0_um * (1 - hlc_band / 2), hlc_lam0_um * (1 + hlc_band / 2), 3)
    hlc_lams_m = hlc_lams_um * 1.0e-6
    spec_lam0_um = 0.73
    spec_band = 0.15
    spec_lams_um = np.linspace(spec_lam0_um * (1 - spec_band / 2), spec_lam0_um * (1 + spec_band / 2), 3)
    spec_lams_m = spec_lams_um * 1.0e-6
    spec_short_lam0_um = 0.66
    spec_short_band = 0.15
    spec_short_lams_um = np.linspace(spec_short_lam0_um * (1 - spec_short_band / 2), spec_short_lam0_um * (1 + spec_short_band / 2), 3)
    spec_short_lams_m = spec_short_lams_um * 1.0e-6
    wide_lam0_um = 0.825
    wide_band = 0.1
    wide_lams_um = np.linspace(wide_lam0_um * (1 - wide_band / 2), wide_lam0_um * (1 + wide_band / 2), 3)
    wide_lams_m = wide_lams_um * 1.0e-6
    return {
        "compact_hlc": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": hlc_lams_m,
            "passvalue": {
                "cor_type": "hlc",
                "use_hlc_dm_patterns": 1,
                "final_sampling_lam0": 0.1,
            },
            "description": "Compact HLC model with default DM patterns over 10% band",
        },
        "full_hlc": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": hlc_lams_m,
            "passvalue": {
                "cor_type": "hlc",
                "use_hlc_dm_patterns": 1,
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
            },
            "description": "Full HLC model with default DM patterns over 10% band",
        },
        "compact_spc_spec_long": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": spec_lams_m,
            "passvalue": {
                "cor_type": "spc-spec_long",
                "final_sampling_lam0": 0.1,
            },
            "description": "Compact SPC spec-long model over 15% band",
        },
        "compact_spc_spec_short": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": spec_short_lams_m,
            "passvalue": {
                "cor_type": "spc-spec_short",
                "final_sampling_lam0": 0.1,
            },
            "description": "Compact SPC spec-short model over 15% band",
        },
        "full_spc_spec_long": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": spec_lams_m,
            "passvalue": {
                "cor_type": "spc-spec_long",
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
            },
            "description": "Full SPC spec-long model over 15% band without error maps",
        },
        "full_spc_spec_short": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": spec_short_lams_m,
            "passvalue": {
                "cor_type": "spc-spec_short",
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
            },
            "description": "Full SPC spec-short model over 15% band without error maps",
        },
        "compact_spc_wide": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": wide_lams_m,
            "passvalue": {
                "cor_type": "spc-wide",
                "final_sampling_lam0": 0.1,
            },
            "description": "Compact SPC wide-field model over 10% band",
        },
        "full_spc_wide": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": wide_lams_m,
            "passvalue": {
                "cor_type": "spc-wide",
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
            },
            "description": "Full SPC wide-field model over 10% band without error maps",
        },
        "compact_hlc_source_offset": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": hlc_lams_m,
            "passvalue": {
                "cor_type": "hlc",
                "use_hlc_dm_patterns": 1,
                "final_sampling_lam0": 0.1,
                "source_x_offset": 3.0,
            },
            "description": "Compact HLC model with nonzero lambda-D source offset",
        },
        "full_hlc_no_field_stop": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": hlc_lams_m,
            "passvalue": {
                "cor_type": "hlc",
                "use_hlc_dm_patterns": 1,
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
                "use_field_stop": 0,
            },
            "description": "Full HLC model without field stop",
        },
        "full_spc_spec_long_no_pupil_mask": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": spec_lams_m,
            "passvalue": {
                "cor_type": "spc-spec_long",
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
                "use_pupil_mask": 0,
            },
            "description": "Full SPC spec-long model without pupil mask",
        },
        "full_none": {
            "func": wfirst_phaseb,
            "output_dim": 128,
            "wavelengths_m": hlc_lams_m,
            "passvalue": {
                "cor_type": "none",
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
                "use_fpm": 0,
            },
            "description": "Full pupil-only model without coronagraph elements",
        },
    }


def run_case(case):
    outputs = []
    samplings = []
    for lam_m in case["wavelengths_m"]:
        field, sampling = case["func"](float(lam_m), int(case["output_dim"]), dict(case["passvalue"]))
        outputs.append(np.asarray(field))
        samplings.append(float(sampling))
    return np.stack(outputs), np.asarray(samplings, dtype=np.float64)


def summarize_output(stack, samplings):
    mags = np.abs(stack)
    return {
        "dtype": str(stack.dtype),
        "shape": list(stack.shape),
        "sampling": [float(x) for x in samplings.tolist()],
        "peak_abs": float(np.max(mags)),
        "mean_abs": float(np.mean(mags)),
        "sum_abs2": float(np.sum(mags * mags)),
    }


def write_output_stack(prefix: str, stack: np.ndarray) -> None:
    prefix_path = Path(prefix)
    prefix_path.parent.mkdir(parents=True, exist_ok=True)
    for idx, plane in enumerate(stack, start=1):
        fits.writeto(f"{prefix}_{idx}_real.fits", np.asarray(np.real(plane), dtype=np.float64), overwrite=True)
        fits.writeto(f"{prefix}_{idx}_imag.fits", np.asarray(np.imag(plane), dtype=np.float64), overwrite=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark external Python WFIRST Phase B models")
    parser.add_argument(
        "--case",
        default="compact_hlc",
    )
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--write-output-prefix", default="")
    parser.add_argument("--parity-only", action="store_true")
    args = parser.parse_args()

    proper, wfirst_phaseb, wfirst_phaseb_compact, proper_root, models_python_root, data_root = load_python_models()
    cases = case_definitions(wfirst_phaseb, wfirst_phaseb_compact)
    if args.case not in cases:
        raise SystemExit(f"unsupported case {args.case}")
    case = cases[args.case]

    if data_root is None:
        write_skipped_report(args.case, "Missing WFIRST Phase B data tree. Set WFIRST_PHASEB_DATA_ROOT to a checkout containing data_phaseb.")
        return

    outdir = project_root() / "bench" / "reports"
    outdir.mkdir(parents=True, exist_ok=True)

    stack, samplings = run_case(case)
    if args.write_output_prefix:
        write_output_stack(args.write_output_prefix, stack)

    timed = (not args.parity_only) and args.samples > 0
    samples_ns = []
    if timed:
        for _ in range(args.samples):
            t0 = time.perf_counter_ns()
            run_case(case)
            samples_ns.append(time.perf_counter_ns() - t0)

    report = {
        "meta": {
            "run_tag": f"wfirst_phaseb_python_{args.case}",
            "backend": "cpu",
            "python_version": os.sys.version,
            "proper_version": getattr(proper, "__version__", "unknown"),
            "proper_root": str(proper_root),
            "models_python_root": str(models_python_root),
            "data_root": str(data_root),
            "baseline": "wfirst_phaseb_python",
            "case": args.case,
            "description": case["description"],
            "output_dim": case["output_dim"],
            "wavelengths_um": [float(v * 1.0e6) for v in case["wavelengths_m"]],
            "available": True,
        },
        "stats": {
            "timed": timed,
            "median_ns": int(statistics.median(samples_ns)) if timed else None,
            "mean_ns": int(statistics.fmean(samples_ns)) if timed else None,
            "min_ns": int(min(samples_ns)) if timed else None,
            "max_ns": int(max(samples_ns)) if timed else None,
            "samples": len(samples_ns),
        },
        "output": summarize_output(stack, samplings),
        "policy": (
            "external WFIRST Phase B Python model timing; TTFx excluded; direct prescription calls with runtime data_dir override"
            if timed
            else "external WFIRST Phase B Python model parity-only run; no timing samples collected"
        ),
    }

    outpath = outdir / f"python_wfirst_phaseb_{args.case}.json"
    outpath.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
