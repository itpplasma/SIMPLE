#!/usr/bin/env python3
"""
Diagnostic plot of the NCSX VMEC equilibrium with coils overlaid.

This script:
  - Ensures the NCSX VMEC equilibrium and coils files exist in a repo-relative
    test-data directory (golden_record/test_data).
  - Downloads and converts the STELLOPT coils file to SIMPLE format if needed.
  - Plots a poloidal cross-section of the VMEC surface together with the
    corresponding coil cross-section.

Output:
  golden_record/test_data/ncsx_equilibrium_with_coils.png
"""

import math
import os
import sys
from typing import Tuple

import numpy as np

try:
  import netCDF4 as nc
except ImportError:
  print("Error: netCDF4 is required to read wout_ncsx.nc", file=sys.stderr)
  sys.exit(1)

try:
  import matplotlib.pyplot as plt
except ImportError:
  print("Error: matplotlib is required to generate plots", file=sys.stderr)
  sys.exit(1)

try:
  from urllib.request import urlretrieve
except ImportError:
  print("Error: urllib.request not available for downloads", file=sys.stderr)
  sys.exit(1)


NCSX_WOUT_URL = (
  "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/"
  "wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc"
)
NCSX_COILS_URL = (
  "https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00"
)


def _project_root() -> str:
  """Return SIMPLE project root (one level above this script)."""
  here = os.path.abspath(os.path.dirname(__file__))
  return os.path.abspath(os.path.join(here, os.pardir))


def ensure_ncsx_files(test_data_dir: str) -> Tuple[str, str]:
  """
  Ensure NCSX VMEC and coils files exist in test_data_dir.

  Returns paths to (wout_ncsx.nc, coils.c09r00.simple).
  """
  os.makedirs(test_data_dir, exist_ok=True)

  wout_path = os.path.join(test_data_dir, "wout_ncsx.nc")
  coils_simple_path = os.path.join(test_data_dir, "coils.c09r00.simple")
  coils_stellopt_path = os.path.join(test_data_dir, "coils.c09r00")

  if not os.path.exists(wout_path):
    print(f"Downloading NCSX wout: {NCSX_WOUT_URL}")
    urlretrieve(NCSX_WOUT_URL, wout_path)

  if not os.path.exists(coils_simple_path):
    if not os.path.exists(coils_stellopt_path):
      print(f"Downloading NCSX coils: {NCSX_COILS_URL}")
      urlretrieve(NCSX_COILS_URL, coils_stellopt_path)
    print("Converting coils.c09r00 -> coils.c09r00.simple")
    convert_coils_to_simple(coils_stellopt_path, coils_simple_path)

  return wout_path, coils_simple_path


def convert_coils_to_simple(src_path: str, dst_path: str) -> None:
  """
  Convert STELLOPT coils file to SIMPLE coils format.

  SIMPLE format:
    line 1: integer N (number of points)
    next N lines: x  y  z  I
  """
  coords = []
  currents = []
  in_filament = False

  with open(src_path, "r") as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      lower = line.lower()
      if "begin filament" in lower:
        in_filament = True
        continue
      if not in_filament or line.startswith("#"):
        continue
      parts = line.split()
      if len(parts) < 4:
        continue
      try:
        x, y, z, I = map(float, parts[:4])
      except ValueError:
        continue
      coords.append((x, y, z))
      currents.append(I)

  if not coords:
    raise RuntimeError(
      f"No coil points found while converting {src_path} to SIMPLE format"
    )

  with open(dst_path, "w") as f:
    f.write(f"{len(coords)}\n")
    for (x, y, z), I in zip(coords, currents):
      f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")


def load_vmec_poloidal_section(wout_path: str, s_index: int = None, n_theta: int = 400):
  """
  Load a poloidal cross-section (R(θ), Z(θ)) at φ = 0 from a VMEC wout file.

  Uses standard Fourier representation:
    R(s, θ, φ) = sum_mn rmnc(s, mn) * cos(m θ - n Nfp φ)
    Z(s, θ, φ) = sum_mn zmns(s, mn) * sin(m θ - n Nfp φ)
  At φ = 0 the dependence on n drops out.
  """
  with nc.Dataset(wout_path, "r") as ds:
    rmnc = ds.variables["rmnc"][:]  # (radius, mn_mode)
    zmns = ds.variables["zmns"][:]  # (radius, mn_mode)
    xm = ds.variables["xm"][:]      # (mn_mode,)
    # xn is not needed at φ = 0

  ns = rmnc.shape[0]
  if s_index is None:
    s_index = ns // 2

  if s_index < 0 or s_index >= ns:
    raise ValueError(f"s_index {s_index} out of range [0, {ns-1}]")

  rmnc_s = rmnc[s_index, :]
  zmns_s = zmns[s_index, :]

  theta = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=True)
  # angles shape: (n_theta, mn_mode)
  angles = np.outer(theta, xm)

  R = np.cos(angles) @ rmnc_s
  Z = np.sin(angles) @ zmns_s

  return R, Z


def load_coils_cross_section(coils_simple_path: str, phi_tolerance: float = 0.05):
  """
  Load coil points from a SIMPLE coils file and return a cross-section near φ = 0.

  Returns arrays (R_coils, Z_coils).
  """
  with open(coils_simple_path, "r") as f:
    header = f.readline()
    try:
      npts = int(header.strip().split()[0])
    except Exception:
      raise RuntimeError(
        f"First line of {coils_simple_path!r} must contain number of points"
      )

    xs = []
    ys = []
    zs = []
    for i in range(npts):
      line = f.readline()
      if not line:
        break
      parts = line.split()
      if len(parts) < 3:
        continue
      try:
        x, y, z = map(float, parts[:3])
      except ValueError:
        continue
      xs.append(x)
      ys.append(y)
      zs.append(z)

  xs = np.asarray(xs, dtype=float)
  ys = np.asarray(ys, dtype=float)
  zs = np.asarray(zs, dtype=float)

  R = np.sqrt(xs * xs + ys * ys)
  phi = np.arctan2(ys, xs)
  Z = zs

  # Select points near φ = 0 (mod 2π)
  phi_mod = np.mod(phi + math.pi, 2.0 * math.pi) - math.pi
  mask = np.abs(phi_mod) < phi_tolerance

  return R[mask], Z[mask]


def main(argv=None):
  if argv is None:
    argv = sys.argv[1:]

  root = _project_root()
  if len(argv) >= 1:
    test_data_dir = os.path.abspath(argv[0])
  else:
    test_data_dir = os.path.join(root, "golden_record", "test_data")

  print(f"Using test data directory: {test_data_dir}")

  wout_path, coils_simple_path = ensure_ncsx_files(test_data_dir)

  print(f"VMEC file:  {wout_path}")
  print(f"Coils file: {coils_simple_path}")

  R_vmec, Z_vmec = load_vmec_poloidal_section(wout_path)
  R_coils, Z_coils = load_coils_cross_section(coils_simple_path)

  fig, ax = plt.subplots(figsize=(6, 6))
  ax.plot(R_vmec, Z_vmec, "b-", label="VMEC surface (s ~ mid)")
  if R_coils.size > 0:
    ax.scatter(R_coils, Z_coils, s=10, c="r", marker="o", label="Coils cross-section")
  else:
    print("Warning: no coil points near φ = 0 found for cross-section")

  ax.set_xlabel("R")
  ax.set_ylabel("Z")
  ax.set_aspect("equal", adjustable="box")
  ax.grid(True, alpha=0.3)
  ax.set_title("NCSX VMEC equilibrium with coils (poloidal cross-section)")
  ax.legend()

  output_path = os.path.join(test_data_dir, "ncsx_equilibrium_with_coils.png")
  plt.tight_layout()
  plt.savefig(output_path, dpi=150)
  plt.close(fig)

  print(f"Diagnostic plot written to: {output_path}")


if __name__ == "__main__":
  main()

