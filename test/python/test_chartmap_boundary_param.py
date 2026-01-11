from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


@pytest.mark.network
def test_chartmap_boundary_param_wall_area(tmp_path: Path) -> None:
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")
    pytest.importorskip("libneo")

    from libneo.chartmap import write_chartmap_from_vmec_boundary
    from pysimple.engineering import compute_wall_area_from_chartmap

    wout_url = "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    theta_chart = tmp_path / "chartmap_theta.nc"
    arc_chart = tmp_path / "chartmap_arc.nc"

    write_chartmap_from_vmec_boundary(
        wout,
        theta_chart,
        nrho=33,
        ntheta=65,
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="theta",
    )
    write_chartmap_from_vmec_boundary(
        wout,
        arc_chart,
        nrho=33,
        ntheta=65,
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="arc",
    )

    area_theta = float(compute_wall_area_from_chartmap(theta_chart))
    area_arc = float(compute_wall_area_from_chartmap(arc_chart))

    assert area_theta > 0.0
    assert area_arc > 0.0

    rel_diff = abs(area_arc - area_theta) / area_theta
    assert rel_diff < 5.0e-3
