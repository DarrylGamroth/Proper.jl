#!/usr/bin/env python3
import argparse
import html
import json
import re
import shutil
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile

ARCHIVE_NAME = "roman_preflight_proper_public_v2.0.2_python.zip"
DOWNLOAD_PAGE = f"https://sourceforge.net/projects/cgisim/files/{ARCHIVE_NAME}/download"
CACHE_DIRNAME = ".cache/wfirst"
OUTPUT_DIRNAME = "data_phaseb_from_roman_preflight"
OLD_LAM_OCCS = [
    "5.4625e-07", "5.49444444444e-07", "5.52638888889e-07", "5.534375e-07", "5.55833333333e-07", "5.59027777778e-07",
    "5.60625e-07", "5.62222222222e-07", "5.65416666667e-07", "5.678125e-07", "5.68611111111e-07", "5.71805555556e-07",
    "5.75e-07", "5.78194444444e-07", "5.81388888889e-07", "5.821875e-07", "5.84583333333e-07", "5.87777777778e-07",
    "5.89375e-07", "5.90972222222e-07", "5.94166666667e-07", "5.965625e-07", "5.97361111111e-07", "6.00555555556e-07", "6.0375e-07",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def cache_root() -> Path:
    return project_root() / CACHE_DIRNAME


def archive_path() -> Path:
    return cache_root() / ARCHIVE_NAME


def compat_root() -> Path:
    return cache_root() / OUTPUT_DIRNAME


def resolve_download_url() -> str:
    request = Request(DOWNLOAD_PAGE, headers={"user-agent": "Mozilla/5.0"})
    with urlopen(request, timeout=30) as response:
        page = response.read().decode("utf-8", errors="replace")
    match = re.search(r"(https://downloads\\.sourceforge\\.net/project/cgisim/" + re.escape(ARCHIVE_NAME) + r"\\?[^\"']+)", page)
    if not match:
        raise RuntimeError("Unable to locate SourceForge mirror URL for Roman preflight archive")
    return html.unescape(match.group(1))


def ensure_archive() -> Path:
    path = archive_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.stat().st_size > 10_000_000:
        return path

    url = resolve_download_url()
    request = Request(url, headers={"user-agent": "Mozilla/5.0"})
    with urlopen(request, timeout=120) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    return path


def extract_bytes(zf: ZipFile, src: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(zf.read(src))


def build_compat_root(force: bool = False) -> Path:
    out = compat_root()
    if out.exists() and not force:
        return out
    if out.exists():
        shutil.rmtree(out)

    zip_path = ensure_archive()
    with ZipFile(zip_path) as zf:
        base = "roman_preflight_proper_public_v2.0.2_python/roman_preflight_proper/preflight_data"
        hlc = f"{base}/hlc_20190210b"
        pol = f"{base}/pol"
        spc_spec = f"{base}/spc_20200617_spec"
        spc_wide = f"{base}/spc_20200610_wfov"

        available_wavelengths = []
        pattern = re.compile(re.escape(hlc) + r"/hlc_fpm_trans_(\d+\.\d+)um_(real|imag)\.fits$")
        for name in zf.namelist():
            match = pattern.match(name)
            if match:
                available_wavelengths.append(float(match.group(1)) * 1.0e-6)
        available_wavelengths = sorted(set(available_wavelengths))
        if not available_wavelengths:
            raise RuntimeError("No HLC FPM transmission files found in Roman preflight archive")

        direct_map = {
            f"{hlc}/pupil_rotated.fits": out / "hlc_20190210" / "run461_pupil_rotated.fits",
            f"{hlc}/lyot.fits": out / "hlc_20190210" / "run461_lyot.fits",
            f"{hlc}/hlc_dm1.fits": out / "hlc_20190210" / "run461_dm1wfe.fits",
            f"{hlc}/hlc_dm2.fits": out / "hlc_20190210" / "run461_dm2wfe.fits",
            f"{hlc}/dm2mask.fits": out / "hlc_20190210" / "run461_dm2mask.fits",
            f"{spc_spec}/pupil_SPC-20200617_1000_rotated.fits": out / "spc_20190130" / "pupil_SPC-20190130_rotated.fits",
            f"{spc_spec}/SPM_SPC-20200617_1000_rounded9.fits": out / "spc_20190130" / "SPM_SPC-20190130.fits",
            f"{spc_spec}/SPM_SPC-20200617_1000_rounded9_rotated.fits": out / "spc_20190130" / "SPM_SPC-20190130_rotated.fits",
            f"{spc_spec}/fpm_0.05lamD.fits": out / "spc_20190130" / "fpm_0.05lamdivD.fits",
            f"{spc_spec}/LS_SPC-20200617_1000.fits": out / "spc_20190130" / "LS_SPC-20190130.fits",
            f"{spc_spec}/LS_SPC-20200617_500.fits": out / "spc_20190130" / "lyotstop_0.5mag.fits",
            f"{spc_wide}/pupil_SPC-20200610_1000_rotated.fits": out / "spc_20181220" / "pupil_SPC-20181220_1k_rotated.fits",
            f"{spc_wide}/SPM_SPC-20200610_1000_rounded9_gray.fits": out / "spc_20181220" / "SPM_SPC-20181220_1000_rounded9_gray.fits",
            f"{spc_wide}/SPM_SPC-20200610_1000_rounded9_gray_rotated.fits": out / "spc_20181220" / "SPM_SPC-20181220_1000_rounded9_gray_rotated.fits",
            f"{spc_wide}/FPM_SPC-20200610_0.1_lamc_div_D.fits": out / "spc_20181220" / "fpm_0.05lamdivD.fits",
            f"{spc_wide}/LS_SPC-20200610_1000.fits": out / "spc_20181220" / "LS_SPC-20181220_1k.fits",
            f"{spc_wide}/LS_SPC-20200610_500.fits": out / "spc_20181220" / "LS_half_symm_CGI180718_Str3.20pct_38D91_N500_pixel.fits",
            f"{pol}/new_toma_amp.fits": out / "pol" / "new_toma_amp.fits",
            f"{pol}/new_toma_pha.fits": out / "pol" / "new_toma_pha.fits",
        }
        for src, dst in direct_map.items():
            extract_bytes(zf, src, dst)

        for lam_label in OLD_LAM_OCCS:
            requested_lam_m = float(lam_label)
            nearest_lam_m = min(available_wavelengths, key=lambda candidate: abs(candidate - requested_lam_m))
            stem = f"hlc_fpm_trans_{nearest_lam_m * 1.0e6:0.8f}um"
            real_bytes = zf.read(f"{hlc}/{stem}_real.fits")
            imag_bytes = zf.read(f"{hlc}/{stem}_imag.fits")
            for suffix, payload in (
                ("real.fits", real_bytes),
                ("imag.fits", imag_bytes),
                ("real_rotated.fits", real_bytes),
                ("imag_rotated.fits", imag_bytes),
            ):
                dst = out / "hlc_20190210" / f"run461_occ_lam{lam_label}theta6.69polp_{suffix}"
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_bytes(payload)

    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a local Phase B compatibility data root from the public Roman preflight archive")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = build_compat_root(force=args.force)
    report = {
        "compat_root": str(root.resolve()),
        "archive": str(archive_path().resolve()),
        "fits_count": sum(1 for _ in root.rglob("*.fits")),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
