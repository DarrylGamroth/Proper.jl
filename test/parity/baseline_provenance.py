#!/usr/bin/env python3
import hashlib
import importlib.metadata
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


SOURCEFORGE_BASELINE_URL = (
    "https://sourceforge.net/projects/proper-library/files/proper_v3.3.4_python.zip"
)
SOURCE_SUFFIXES = {".c", ".fits", ".h", ".py"}


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_snapshot_sha256(pyproper_root):
    root = Path(pyproper_root).resolve()
    digest = hashlib.sha256()
    paths = sorted(
        path
        for path in (root / "proper").rglob("*")
        if path.is_file() and path.suffix.lower() in SOURCE_SUFFIXES
    )
    for path in paths:
        relative = path.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(_sha256_file(path)))
    return digest.hexdigest()


def generated_at_utc():
    source_date_epoch = os.environ.get("SOURCE_DATE_EPOCH")
    instant = (
        datetime.fromtimestamp(int(source_date_epoch), tz=timezone.utc)
        if source_date_epoch is not None
        else datetime.now(tz=timezone.utc)
    )
    return instant.isoformat().replace("+00:00", "Z")


def installed_version(package):
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def write_provenance(
    output_path,
    *,
    generator,
    pyproper_root,
    proper,
    seed,
    cases,
    artifacts,
):
    output_path = Path(output_path)
    artifact_root = output_path.parent
    artifact_hashes = {
        name: _sha256_file(artifact_root / name) for name in sorted(artifacts)
    }
    payload = {
        "schema_version": 1,
        "baseline": {
            "name": "PROPER Python",
            "version": str(getattr(proper, "__version__", "3.3.4")),
            "source_url": SOURCEFORGE_BASELINE_URL,
            "source_root": str(Path(pyproper_root).resolve()),
            "source_snapshot_sha256": source_snapshot_sha256(pyproper_root),
        },
        "generator": generator,
        "generated_at_utc": generated_at_utc(),
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "numpy": installed_version("numpy"),
            "scipy": installed_version("scipy"),
            "astropy": installed_version("astropy"),
            "backend": "NumPy/SciPy CPU",
        },
        "numeric_precision": ["float64", "complex128"],
        "rng_seed": seed,
        "cases": cases,
        "artifacts": artifact_hashes,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
