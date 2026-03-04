#!/usr/bin/env python3
import json
import os
import statistics
import time
from pathlib import Path

import numpy as np


def workload(n=512):
    a = np.ones((n, n), dtype=np.complex128)
    x = np.linspace(-1.0, 1.0, n)
    xx, yy = np.meshgrid(x, x)
    phase = np.exp(1j * np.pi * (xx * xx + yy * yy))
    a *= phase
    psf = np.abs(a) ** 2
    return psf


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "reports"
    outdir.mkdir(parents=True, exist_ok=True)

    # Warmup for fair steady-state timings.
    workload()

    samples_ns = []
    for _ in range(50):
        t0 = time.perf_counter_ns()
        workload()
        samples_ns.append(time.perf_counter_ns() - t0)

    report = {
        "meta": {
            "run_tag": "steady_state",
            "backend": "cpu",
            "python_version": os.sys.version,
        },
        "stats": {
            "median_ns": int(statistics.median(samples_ns)),
            "mean_ns": int(statistics.fmean(samples_ns)),
            "min_ns": int(min(samples_ns)),
            "max_ns": int(max(samples_ns)),
            "samples": len(samples_ns),
        },
        "policy": "steady-state timing only",
    }

    (outdir / "python_steady_state.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
