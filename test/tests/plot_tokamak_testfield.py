#!/usr/bin/env python3
"""Plot the canonical tokamak TEST-field orbit used in symplectic benchmarks."""

from pathlib import Path
import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as exc:
    # Allow this script to be used as an optional plot helper in test
    # environments where matplotlib/numpy may not be installed.
    print(f"Skipping tokamak TEST-field plot (missing dependency: {exc})")
    raise SystemExit(0)


OUT_FILE = Path("tokamak_testfield_orbit.png")
def load_orbit(path: Path):
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError(f"Unexpected layout in {path}")
    r = data[:, 0]
    theta = data[:, 1]
    h = data[:, 4]
    R = 1.0 + r * np.cos(theta)
    Z = r * np.sin(theta)
    return R, Z, h


def main() -> int:
    files = sorted(Path(".").glob("tokamak_testfield_*.dat"))
    if not files:
        print("No tokamak_testfield_*.dat files found; run test_sympl_tok.x first.")
        return 1

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].set_xlabel("R/R0")
    axes[0].set_ylabel("Z/R0")
    axes[0].set_title("Tokamak TEST field orbits")

    axes[1].set_xlabel("Step")
    axes[1].set_ylabel("H/H0")
    axes[1].set_title("Hamiltonian drift")

    markers = ["o", "s", "^", "x", "d", "v"]

    for idx, path in enumerate(files):
        R, Z, h = load_orbit(path)
        label = path.stem.replace("tokamak_testfield_", "")
        marker = markers[idx % len(markers)]
        axes[0].plot(R, Z, marker=marker, linestyle="none", markersize=2, label=label)
        axes[1].plot(
            h / h[0],
            marker=marker,
            linestyle="none",
            markersize=2,
            label=label,
        )

    axes[0].legend()
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(OUT_FILE, dpi=150)
    print(f"Saved {OUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
