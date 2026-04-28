#!/usr/bin/env python3
import json
import os
import statistics
import sys
import time
from pathlib import Path

GRID_N = 512
N_SAMPLES = 20


def load_proper():
    project_root = Path(__file__).resolve().parents[2]
    proper_root = Path(os.environ.get("PYPROPER_ROOT", project_root / ".." / "proper_v3.3.4_python")).resolve()
    if not proper_root.is_dir():
        raise RuntimeError(f"Missing Python baseline source tree: {proper_root}")

    sys.path.insert(0, str(proper_root))
    try:
        import proper  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - benchmark harness guardrail
        raise RuntimeError(
            "Unable to import Python PROPER baseline. "
            f"Interpreter: {sys.executable}. "
            "Use .venv-parity with numpy/scipy/astropy/matplotlib installed "
            "or run ./scripts/setup_parity_venv.sh."
        ) from exc

    proper.print_it = False
    proper.verbose = False
    proper.print_total_intensity = False
    proper.do_table = False

    return proper, proper_root


def workload(proper, n=GRID_N):
    wf = proper.prop_begin(2.4, 0.55e-6, n, beam_diam_fraction=0.5)
    proper.prop_circular_aperture(wf, 0.6)
    proper.prop_lens(wf, 20.0)
    proper.prop_propagate(wf, 20.0)
    proper.prop_end(wf)


def main() -> None:
    proper, proper_root = load_proper()
    outdir = Path(__file__).resolve().parents[1] / "reports"
    outdir.mkdir(parents=True, exist_ok=True)

    workload(proper)

    samples_ns = []
    for _ in range(N_SAMPLES):
        t0 = time.perf_counter_ns()
        workload(proper)
        samples_ns.append(time.perf_counter_ns() - t0)

    report = {
        "meta": {
            "run_tag": "steady_state",
            "backend": "cpu",
            "python_version": os.sys.version,
            "proper_version": getattr(proper, "__version__", "unknown"),
            "proper_root": str(proper_root),
            "grid_n": GRID_N,
            "baseline": "python334_patched",
        },
        "stats": {
            "median_ns": int(statistics.median(samples_ns)),
            "mean_ns": int(statistics.fmean(samples_ns)),
            "min_ns": int(min(samples_ns)),
            "max_ns": int(max(samples_ns)),
            "samples": len(samples_ns),
        },
        "policy": "steady-state PROPER workload timing only; TTFx excluded",
    }

    (outdir / "python_steady_state.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
