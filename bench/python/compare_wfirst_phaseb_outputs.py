#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def loadjson(path: Path):
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def timed_median_ms(report):
    if report is None:
        return None
    stats = report.get("stats", {})
    if not stats.get("timed", True):
        return None
    median_ns = stats.get("median_ns")
    if median_ns is None:
        return None
    return float(median_ns) / 1.0e6


def ratio_or_none(a, b):
    if a is None or b is None:
        return None
    return a / b


def load_case_stack(prefix: Path):
    planes = []
    idx = 1
    while (prefix.parent / f"{prefix.name}_{idx}_real.fits").is_file():
        real = fits.getdata(prefix.parent / f"{prefix.name}_{idx}_real.fits")
        imag = fits.getdata(prefix.parent / f"{prefix.name}_{idx}_imag.fits")
        planes.append(real + 1j * imag)
        idx += 1
    if not planes:
        raise RuntimeError(f"no FITS outputs found for prefix {prefix}")
    return np.stack(planes, axis=0)


def case_comparison(case_name: str):
    root = project_root() / "bench" / "reports"
    py_report = loadjson(root / f"python_wfirst_phaseb_{case_name}.json")
    jl_report = loadjson(root / f"julia_wfirst_phaseb_{case_name}.json")
    if py_report is None or jl_report is None:
        return {
            "case": case_name,
            "available": False,
            "python_median_ms": None,
            "julia_median_ms": None,
            "python_over_julia": None,
            "relative_l2": None,
            "max_abs_diff": None,
            "max_sampling_abs_diff": None,
            "threads": None,
            "python_timed": False,
            "julia_timed": False,
            "reason": "missing report artifact",
        }

    py_prefix = root / f"python_wfirst_phaseb_{case_name}"
    jl_prefix = root / f"julia_wfirst_phaseb_{case_name}"
    if not (root / f"python_wfirst_phaseb_{case_name}_1_real.fits").is_file() or not (root / f"julia_wfirst_phaseb_{case_name}_1_real.fits").is_file():
        return {
            "case": case_name,
            "available": False,
            "python_median_ms": timed_median_ms(py_report),
            "julia_median_ms": timed_median_ms(jl_report),
            "python_over_julia": ratio_or_none(timed_median_ms(py_report), timed_median_ms(jl_report)),
            "relative_l2": None,
            "max_abs_diff": None,
            "max_sampling_abs_diff": None,
            "threads": int(jl_report["meta"]["threads"]),
            "python_timed": bool(py_report.get("stats", {}).get("timed", True)),
            "julia_timed": bool(jl_report.get("stats", {}).get("timed", True)),
            "reason": "missing FITS output artifact",
        }

    py_stack = load_case_stack(py_prefix)
    jl_stack = load_case_stack(jl_prefix)
    diff = jl_stack - py_stack
    py_norm = np.linalg.norm(py_stack.ravel())
    rel_l2 = np.linalg.norm(diff.ravel()) if py_norm == 0 else np.linalg.norm(diff.ravel()) / py_norm
    max_abs = np.max(np.abs(diff))

    py_sampling = np.asarray(py_report["output"]["sampling"], dtype=np.float64)
    jl_sampling = np.asarray(jl_report["output"]["sampling"], dtype=np.float64)
    sampling_abs = np.max(np.abs(jl_sampling - py_sampling))
    py_ms = timed_median_ms(py_report)
    jl_ms = timed_median_ms(jl_report)

    return {
        "case": case_name,
        "available": True,
        "python_median_ms": py_ms,
        "julia_median_ms": jl_ms,
        "python_over_julia": ratio_or_none(py_ms, jl_ms),
        "relative_l2": float(rel_l2),
        "max_abs_diff": float(max_abs),
        "max_sampling_abs_diff": float(sampling_abs),
        "threads": int(jl_report["meta"]["threads"]),
        "python_timed": bool(py_report.get("stats", {}).get("timed", True)),
        "julia_timed": bool(jl_report.get("stats", {}).get("timed", True)),
    }


def fmt_ms(x):
    return "n/a" if x is None else f"{x:.2f} ms"


def fmt_ratio(x):
    return "n/a" if x is None else f"{x:.2f}x"


def fmt_sci(x):
    return "n/a" if x is None else f"{x:.4g}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="compact_hlc,full_hlc,compact_spc_spec_short,full_spc_spec_short,compact_spc_ifs_short,full_spc_ifs_short,compact_spc_spec_long,full_spc_spec_long,compact_spc_ifs_long,full_spc_ifs_long,compact_spc_wide,full_spc_wide")
    args = parser.parse_args()

    cases = [case for case in args.cases.split(",") if case]
    rows = [case_comparison(case) for case in cases]
    outpath = project_root() / "bench" / "reports" / "wfirst_phaseb_cpu_comparison.json"
    outpath.write_text(json.dumps({"cases": rows}))

    print("WFIRST Phase B CPU Comparison")
    print("============================")
    print(f"{'Case':16}{'Python':>12}{'Julia':>12}{'Py/Jl':>10}{'RelL2':>14}{'MaxAbs':>14}{'dSampling':>14}{'Threads':>10}")
    for row in rows:
        line = (
            f"{row['case']:16}"
            f"{fmt_ms(row['python_median_ms']):>12}"
            f"{fmt_ms(row['julia_median_ms']):>12}"
            f"{fmt_ratio(row['python_over_julia']):>10}"
            f"{fmt_sci(row['relative_l2']):>14}"
            f"{fmt_sci(row['max_abs_diff']):>14}"
            f"{fmt_sci(row['max_sampling_abs_diff']):>14}"
            f"{str(row['threads']) if row['threads'] is not None else 'n/a':>10}"
        )
        if row.get("reason") is not None:
            line += f"  # {row['reason']}"
        print(line)


if __name__ == "__main__":
    main()
