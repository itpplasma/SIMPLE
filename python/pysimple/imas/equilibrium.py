"""Axisymmetric equilibrium IDS to SIMPLE's libneo Boozer chartmap input."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

from ._common import array, check_ids, ids_version, scalar, select_slice


def _geqdsk_input(ids, time=None):
    """Assemble the COCOS-3 intermediary expected by libneo's field reader."""
    from scipy.interpolate import CubicSpline

    check_ids(ids, "equilibrium")
    eq = select_slice(ids.time_slice, time)
    rectangular = [p for p in eq.profiles_2d if int(p.grid_type.index) == 1]
    if len(rectangular) != 1:
        raise ValueError("Provide exactly one rectangular (R,Z) equilibrium grid")
    p2 = rectangular[0]
    r = array(p2.grid.dim1, "equilibrium R grid")
    z = array(p2.grid.dim2, "equilibrium Z grid")
    for grid, name in [(r, "R"), (z, "Z")]:
        if len(grid) < 6 or np.any(np.diff(grid) <= 0):
            raise ValueError(f"{name} grid must have >=6 increasing points")
        if not np.allclose(np.diff(grid), np.diff(grid)[0], rtol=1e-8, atol=1e-12):
            raise ValueError(f"{name} grid must be uniform for GEQDSK conversion")
    if np.any(r <= 0):
        raise ValueError("R must be positive")
    psi = array(p2.psi, "equilibrium psi", ndim=2)
    if psi.shape != (len(r), len(z)):
        raise ValueError(
            "equilibrium psi must have shape (nR, nZ), without transposition"
        )
    p1 = eq.profiles_1d
    psi1 = array(p1.psi, "equilibrium profiles_1d.psi")
    if len(psi1) < 4:
        raise ValueError("Provide >=4 equilibrium profile points")
    delta = np.diff(psi1)
    if not (np.all(delta > 0) or np.all(delta < 0)):
        raise ValueError("equilibrium profile psi must be strictly monotonic")
    order = np.argsort(psi1)
    psi1 = psi1[order]
    q = eq.global_quantities
    axis = scalar(q.psi_axis, "psi_axis")
    edge = scalar(q.psi_boundary, "psi_boundary")
    if axis == edge or min(axis, edge) < psi1[0] or max(axis, edge) > psi1[-1]:
        raise ValueError(
            "Equilibrium profiles must cover distinct axis/boundary fluxes"
        )
    new_psi = np.linspace(axis, edge, len(r))

    def interpolate(node, name):
        values = array(node, name)
        if len(values) != len(psi1):
            raise ValueError(f"{name} length differs from the profile psi grid")
        return CubicSpline(psi1, values[order])

    f = interpolate(p1.f, "equilibrium f")
    pressure = interpolate(p1.pressure, "equilibrium pressure")
    safety = interpolate(p1.q, "equilibrium q")
    # DD3=COCOS11 (+Bpol); DD4=COCOS17 (-Bpol). libneo expects COCOS3:
    # B_R=-psi_Z/R, B_Z=psi_R/R and flux divided by 2*pi, with opposite q.
    scale = (-1 if ids_version(ids).startswith("3.") else 1) / (2 * np.pi)
    outline = eq.boundary.outline
    boundary = np.column_stack(
        (array(outline.r, "boundary R"), array(outline.z, "boundary Z"))
    )
    if len(boundary) < 4:
        raise ValueError("Provide a closed equilibrium boundary outline")
    if not np.allclose(boundary[0], boundary[-1], rtol=0, atol=1e-10):
        raise ValueError("Equilibrium boundary outline must be closed")
    r0 = scalar(ids.vacuum_toroidal_field.r0, "vacuum field reference R")
    if r0 <= 0:
        raise ValueError("vacuum field reference R must be positive")
    # F at the magnetic axis is the relevant toroidal field normalization;
    # B0 is evaluated from the supplied F, avoiding inconsistent B0 metadata.
    gfile = dict(
        header="SIMPLE IMAS equilibrium (COCOS 3)",
        nrgr=len(r),
        nzgr=len(z),
        rboxlength=r[-1] - r[0],
        zboxlength=z[-1] - z[0],
        R0=r0,
        rboxleft=r[0],
        zboxmid=(z[-1] + z[0]) / 2,
        Rpsi0=scalar(q.magnetic_axis.r, "magnetic axis R"),
        Zpsi0=scalar(q.magnetic_axis.z, "magnetic axis Z"),
        PsiaxisVs=scale * axis,
        PsiedgeVs=scale * edge,
        Btor_at_R0=float(f(axis)) / r0,
        Ip=scalar(q.ip, "plasma current"),
        fprof=f(new_psi),
        ptotprof=pressure(new_psi),
        fdfdpsiprof=f(new_psi) * f(new_psi, 1) / scale,
        dpressdpsiprof=pressure(new_psi, 1) / scale,
        PsiVs=scale * psi.T,
        qprof=-safety(new_psi),
        npbound=len(boundary),
        nplimiter=0,
        Lcfs=boundary,
        Limiter=np.empty((0, 2)),
    )
    return gfile


def equilibrium_to_chartmap(ids, output, *, time=None, timeout=1200, **resolution):
    """Convert an equilibrium IDS to a physical Boozer chartmap.

    Requires libneo with its efit_to_boozer extension. The converter runs in
    a child process so a Fortran STOP or global converter state cannot corrupt
    the caller. DD 3 and DD 4 flux/sign conventions are mapped explicitly.
    ``resolution`` accepts libneo's nrho, ntheta, nstep, nlabel, ntheta_int,
    nsurfmax and nsurf options; inspect convergence for production studies.
    """
    from libneo.eqdsk_base import write_eqdsk
    import netCDF4

    allowed = {"nrho", "ntheta", "nstep", "nlabel", "ntheta_int", "nsurfmax", "nsurf"}
    if set(resolution) - allowed:
        raise ValueError(
            f"Unknown chartmap options: {sorted(set(resolution) - allowed)}"
        )
    for key, value in resolution.items():
        if not isinstance(value, int) or value < 8:
            raise ValueError(f"{key} must be an integer >=8")
    data = _geqdsk_input(ids, time)
    output = Path(output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    options = dict(resolution)
    options["psimax"] = abs(data["PsiedgeVs"] - data["PsiaxisVs"]) * 1e8
    with tempfile.TemporaryDirectory(prefix="simple-imas-equilibrium-") as tmp:
        gfile = Path(tmp) / "equilibrium.geqdsk"
        candidate = Path(tmp) / "chartmap.nc"
        write_eqdsk(gfile, data)
        command = (
            "import json,sys; "
            "from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap; "
            "convert_eqdsk_to_chartmap(sys.argv[1],sys.argv[2],**json.loads(sys.argv[3]))"
        )
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    command,
                    str(gfile),
                    str(candidate),
                    json.dumps(options),
                ],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("libneo equilibrium conversion timed out") from exc
        if result.returncode or not candidate.is_file():
            raise RuntimeError(
                "libneo equilibrium conversion failed:\n"
                + (result.stdout + result.stderr)[-4000:]
            )
        with netCDF4.Dataset(candidate, "a") as ds:
            ds.imas_dd_version = ids_version(ids)
            ds.imas_cocos = 11 if ids_version(ids).startswith("3.") else 17
            ds.imas_time = scalar(
                select_slice(ids.time_slice, time).time, "equilibrium time"
            )
            ds.eqdsk2chartmap_source = (
                "equilibrium IDS (see imas_dd_version and imas_time)"
            )
        output.write_bytes(candidate.read_bytes())
    return output
