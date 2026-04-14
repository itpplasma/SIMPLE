#!/usr/bin/env python3
"""Generate the QA GVEC chartmap on the fly and run the roundtrip check."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

from boozer_chartmap_artifacts import download_if_missing
from test_e2e_boozer_chartmap import (
    build_gvec_chartmap,
    assert_boozer_chartmap_file,
    run_roundtrip,
)


TEST_DATA = Path(__file__).resolve().parent.parent / "test_data"
QA_WOUT = TEST_DATA / "wout.nc"
QA_WOUT_URL = (
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/"
    "tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)

GVEC_QA_CACHE_ROOT = Path(
    os.environ.get(
        "SIMPLE_GVEC_QA_CACHE_ROOT",
        str(Path.home() / "data" / "SIMPLE" / "gvec_qa_roundtrip"),
    )
).expanduser()


def _cache_signature() -> str:
    digest = hashlib.sha256()
    digest.update(QA_WOUT_URL.encode("utf-8"))
    digest.update(Path(__file__).read_bytes())
    digest.update((Path(__file__).resolve().parents[2] / "tools" / "gvec_to_boozer_chartmap.py").read_bytes())
    return digest.hexdigest()[:16]


def _cache_dir() -> Path:
    return GVEC_QA_CACHE_ROOT / _cache_signature()


def main() -> None:
    workdir = _cache_dir()
    workdir.mkdir(parents=True, exist_ok=True)

    wout = download_if_missing(QA_WOUT, QA_WOUT_URL)
    chartmap = workdir / "boozer_chartmap.nc"
    if chartmap.exists():
        print(f"Reusing cached QA GVEC chartmap from {chartmap}")
    else:
        print(f"Building QA GVEC chartmap in {workdir}")
        chartmap = build_gvec_chartmap(
            wout,
            workdir,
            minimize_tol=1.0e-10,
            max_iter=4,
            total_iter=8,
        )
    assert_boozer_chartmap_file(chartmap)
    run_roundtrip(
        wout,
        chartmap,
        "external",
        workdir / "gvec",
        "QA: VMEC direct vs GVEC chartmap",
        "GVEC chartmap",
        "2.5e-4",
        "1.0e-6",
    )
    print("QA GVEC roundtrip PASS")


if __name__ == "__main__":
    main()
