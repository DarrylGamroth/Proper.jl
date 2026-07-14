#!/usr/bin/env python3
"""Fail if the Python reference can silently select SciPy interpolation."""

import os
from pathlib import Path

import numpy as np

from python334_runtime import load_python334_runtime


def main():
    project = Path(__file__).resolve().parents[2]
    pyproper_root = Path(
        os.environ.get("PYPROPER_ROOT", project / ".." / "proper_v3.3.4_python")
    ).resolve()
    proper, native_runtime = load_python334_runtime(pyproper_root)
    native_runtime.assert_active()

    image = np.arange(64 * 64, dtype=np.float64).reshape(64, 64)
    coordinates = np.linspace(4.125, 58.875, 64, dtype=np.float64)
    serial = proper.prop_cubic_conv(
        image, coordinates, coordinates, THREADED=False, GRID=True
    )
    threaded = proper.prop_cubic_conv(
        image, coordinates, coordinates, THREADED=True, GRID=True
    )
    np.testing.assert_array_equal(threaded, serial)

    xgrid = np.tile(coordinates.reshape(1, -1), (coordinates.size, 1))
    xgrid = xgrid + np.linspace(0.0, 0.3, coordinates.size).reshape(-1, 1)
    ygrid = np.tile(coordinates.reshape(-1, 1), (1, coordinates.size))
    ygrid = ygrid - np.linspace(0.0, 0.2, coordinates.size).reshape(1, -1)
    corrected = proper.prop_cubic_conv(image, xgrid, ygrid, GRID=False)
    original = proper.prop_cubic_conv._proper_native_original
    corrected_reference = original(
        image,
        xgrid,
        ygrid,
        THREADED=False,
        GRID=False,
    )
    np.testing.assert_array_equal(corrected, corrected_reference)

    def scipy_fallback_was_called(*_args, **_kwargs):
        raise AssertionError("SciPy interpolation fallback was called")

    magnify_globals = proper.prop_magnify.__globals__
    previous_fallback = magnify_globals.get("map_coordinates")
    magnify_globals["map_coordinates"] = scipy_fallback_was_called
    try:
        quick = proper.prop_magnify(image, 0.75, 40, QUICK=True)
        sinc = proper.prop_magnify(image, 0.75, 40)
    finally:
        if previous_fallback is None:
            magnify_globals.pop("map_coordinates", None)
        else:
            magnify_globals["map_coordinates"] = previous_fallback

    if quick.shape != (40, 40) or sinc.shape != (40, 40):
        raise AssertionError("Native interpolation smoke test returned the wrong shape")
    if not np.isfinite(quick).all() or not np.isfinite(sinc).all():
        raise AssertionError("Native interpolation smoke test returned non-finite values")
    native_runtime.assert_active()
    print("Pinned Python 3.3.4 native interpolation runtime verified.")


if __name__ == "__main__":
    main()
