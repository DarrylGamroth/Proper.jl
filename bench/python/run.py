#!/usr/bin/env python3
import cmath
import json
import os
import statistics
import time
from pathlib import Path


def workload(n=256):
    # Pure-stdlib CPU workload to avoid optional dependencies in CI/dev environments.
    total = 0.0
    half = n / 2.0
    for i in range(n):
        y = (i - half) / half
        for j in range(n):
            x = (j - half) / half
            phase = cmath.exp(1j * 3.141592653589793 * (x * x + y * y))
            total += (phase.real * phase.real + phase.imag * phase.imag)
    return total


def main() -> None:
    outdir = Path(__file__).resolve().parents[1] / "reports"
    outdir.mkdir(parents=True, exist_ok=True)

    workload()

    samples_ns = []
    for _ in range(20):
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
