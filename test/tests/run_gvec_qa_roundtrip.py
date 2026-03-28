#!/usr/bin/env python3
"""Generate the QA GVEC chartmap on the fly and run the roundtrip check."""

from __future__ import annotations

from pathlib import Path
import shutil

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


def main() -> None:
    workdir = Path.cwd() / "boozer_chartmap_gvec_qa"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    wout = download_if_missing(QA_WOUT, QA_WOUT_URL)
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
