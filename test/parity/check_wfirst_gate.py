#!/usr/bin/env python3
"""Contract checks for the executable WFIRST parity threshold gate."""

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "bench" / "python"))

from compare_wfirst_phaseb_outputs import apply_parity_thresholds  # noqa: E402


def row(**overrides):
    value = {
        "case": "synthetic",
        "available": True,
        "relative_l2": 1e-15,
        "max_abs_diff": 1e-16,
        "max_sampling_abs_diff": 0.0,
    }
    value.update(overrides)
    return value


def main():
    passing = apply_parity_thresholds(row())
    if not passing["pass"] or passing["failures"]:
        raise AssertionError("machine-precision WFIRST comparison did not pass")

    inaccurate = apply_parity_thresholds(row(relative_l2=4.13e-3))
    if inaccurate["pass"] or not inaccurate["failures"]:
        raise AssertionError("historical DM mismatch was not rejected")

    unavailable = apply_parity_thresholds(
        row(available=False, reason="missing fixture")
    )
    if unavailable["pass"] or "missing fixture" not in unavailable["failures"][0]:
        raise AssertionError("unavailable WFIRST case was not rejected")

    nonfinite = apply_parity_thresholds(row(max_abs_diff=float("nan")))
    if nonfinite["pass"] or not nonfinite["failures"]:
        raise AssertionError("non-finite WFIRST metric was not rejected")

    print("WFIRST numerical parity gate contract verified.")


if __name__ == "__main__":
    main()
