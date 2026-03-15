#!/usr/bin/env python3
import argparse
import json
import os
import statistics
import sys
import time
from pathlib import Path

import numpy as np


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
    lam0_um = 0.575
    band = 0.1
    lams_um = np.linspace(lam0_um * (1 - band / 2), lam0_um * (1 + band / 2), 3)
    lams_m = lams_um * 1.0e-6
    return {
        "compact_hlc": {
            "func": wfirst_phaseb_compact,
            "output_dim": 128,
            "wavelengths_m": lams_m,
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
            "wavelengths_m": lams_m,
            "passvalue": {
                "cor_type": "hlc",
                "use_hlc_dm_patterns": 1,
                "final_sampling_lam0": 0.1,
                "use_errors": 0,
            },
            "description": "Full HLC model with default DM patterns over 10% band",
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark external Python WFIRST Phase B models")
    parser.add_argument("--case", default="compact_hlc", choices=("compact_hlc", "full_hlc"))
    parser.add_argument("--samples", type=int, default=3)
    args = parser.parse_args()

    proper, wfirst_phaseb, wfirst_phaseb_compact, proper_root, models_python_root, data_root = load_python_models()
    cases = case_definitions(wfirst_phaseb, wfirst_phaseb_compact)
    case = cases[args.case]

    if data_root is None:
        write_skipped_report(args.case, "Missing WFIRST Phase B data tree. Set WFIRST_PHASEB_DATA_ROOT to a checkout containing data_phaseb.")
        return

    outdir = project_root() / "bench" / "reports"
    outdir.mkdir(parents=True, exist_ok=True)

    stack, samplings = run_case(case)

    samples_ns = []
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
            "median_ns": int(statistics.median(samples_ns)),
            "mean_ns": int(statistics.fmean(samples_ns)),
            "min_ns": int(min(samples_ns)),
            "max_ns": int(max(samples_ns)),
            "samples": len(samples_ns),
        },
        "output": summarize_output(stack, samplings),
        "policy": "external WFIRST Phase B Python model timing only; TTFx excluded; direct prescription calls with runtime data_dir override",
    }

    outpath = outdir / f"python_wfirst_phaseb_{args.case}.json"
    outpath.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
