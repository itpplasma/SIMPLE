"""Read-only plots of SIMPLE output; no command starts an orbit trace."""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import NormalDist

import numpy as np


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "install plotting dependencies with pysimple[plotting]"
        ) from exc
    return plt


def wilson_interval(successes, trials, confidence=0.95):
    """Wilson limits for binomial counts, with NaN limits for zero trials."""
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between 0 and 1")
    k, n = np.broadcast_arrays(
        np.asarray(successes, dtype=float), np.asarray(trials, dtype=float)
    )
    if (
        not np.all(np.isfinite(k))
        or not np.all(np.isfinite(n))
        or np.any(n < 0)
        or np.any(k < 0)
        or np.any(k > n)
    ):
        raise ValueError("counts must be finite and satisfy 0 <= successes <= trials")
    z = NormalDist().inv_cdf(0.5 + confidence / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n
        denominator = 1 + z**2 / n
        centre = (p + z**2 / (2 * n)) / denominator
        half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denominator
    return np.maximum(0, centre - half), np.minimum(1, centre + half)


def _table(path, columns):
    if not path.is_file():
        raise FileNotFoundError(f"required output file not found: {path}")
    try:
        data = np.loadtxt(path, ndmin=2)
    except ValueError as exc:
        raise ValueError(f"could not parse {path}: {exc}") from exc
    if not data.size or data.shape[1] < columns:
        raise ValueError(f"{path} must contain at least {columns} columns and one row")
    return data


def _finish(fig, output, show):
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    display = output is None if show is None else show
    if display:
        _pyplot().show()
    return fig


def plot_losses(directory, output=None, confidence=0.95, show=None):
    """Confinement split and Wilson band for physical losses among resolved markers."""
    directory = Path(directory)
    path = directory / "confined_fraction.dat"
    data = _table(path, 3)
    if not np.all(np.isfinite(data[:, 0])) or np.any(data[:, 0] < 0):
        raise ValueError(f"{path}: times must be finite and nonnegative")
    # t=0 cannot be displayed on a logarithmic axis.
    data = data[data[:, 0] > 0]
    if not len(data):
        raise ValueError(f"{path} contains no positive time samples")
    time, passing, trapped = data[:, :3].T
    if data.shape[1] >= 4:
        trials = data[:, 3]
    else:
        trials = np.full(len(data), len(_table(directory / "times_lost.dat", 2)))
    if (
        not np.all(np.isfinite(time))
        or np.any(np.diff(time) <= 0)
        or not np.all(np.isfinite(trials))
        or np.any(trials < 0)
        or np.any(trials != np.rint(trials))
    ):
        raise ValueError(
            f"{path}: times must increase and marker counts be nonnegative integers"
        )
    valid = trials > 0
    if not np.any(valid):
        raise ValueError(f"{path} contains no resolved markers")
    confined = passing + trapped
    if (
        not np.all(np.isfinite(data[valid, 1:3]))
        or np.any(data[valid, 1:3] < 0)
        or np.any(confined[valid] > 1 + 1e-12)
    ):
        raise ValueError(
            f"{path}: passing/trapped fractions must be finite and sum to at most 1"
        )
    # Preserve gaps for rows with no resolved markers; these are not physical losses.
    passing = np.where(valid, passing, np.nan)
    trapped = np.where(valid, trapped, np.nan)
    confined = np.where(valid, np.minimum(confined, 1), np.nan)
    lost = 1 - confined
    counts = np.rint(np.where(valid, lost, 0) * trials)
    lower, upper = wilson_interval(counts, trials, confidence)
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.semilogx(time, confined, color="black", label="Confined")
    # Matplotlib's blue/orange pair, also distinguished by line style.
    ax.semilogx(time, passing, color="tab:blue", linestyle="-.", label="Passing")
    ax.semilogx(time, trapped, color="tab:orange", linestyle=":", label="Trapped")
    ax.semilogx(time, lost, color="tab:red", linestyle="--", label="Lost")
    ax.fill_between(
        time,
        lower,
        upper,
        color="tab:red",
        alpha=0.22,
        label=f"Lost {confidence:.0%} Wilson interval",
    )
    ax.set(
        xlabel="Time / s",
        ylabel="Fraction of resolved markers",
        ylim=(0, 1.05),
        title="SIMPLE confinement and loss",
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    return _finish(fig, output, show)


def _results_path(value):
    path = Path(value)
    return path / "results.nc" if path.is_dir() else path


def _netcdf(path, names):
    """Read only requested variables, using named dimensions to order particle rows."""
    if not path.is_file():
        raise FileNotFoundError(f"required NetCDF results file not found: {path}")
    try:
        from netCDF4 import Dataset
    except ImportError as exc:
        raise ImportError("NetCDF plots require pysimple[plotting]") from exc
    values = {}
    with Dataset(path, "r") as ds:
        for name in names:
            if name not in ds.variables:
                raise ValueError(f"{path} has no {name} data")
            var = ds.variables[name]
            data = np.ma.asarray(var[:], dtype=float).filled(np.nan)
            if "particle" not in var.dimensions:
                raise ValueError(f"{path}: {name} has no particle dimension")
            values[name] = np.moveaxis(data, var.dimensions.index("particle"), 0)
        units = getattr(ds, "cart_units", "Cartesian units")
    return values, units


def plot_classes(value, output=None, show=None):
    """Counts of the J_parallel, topology and fractal class codes (0=unclassified)."""
    path = Path(value)
    if path.is_dir() and (path / "class_parts.dat").is_file():
        path = path / "class_parts.dat"
        classes = _table(path, 6)[:, 3:6]
    else:
        path = _results_path(path)
        variables, _ = _netcdf(path, ("iclass",))
        classes = variables["iclass"]
    if (
        classes.ndim != 2
        or classes.shape[1] != 3
        or not classes.size
        or not np.all(np.isfinite(classes))
        or np.any(classes != np.rint(classes))
    ):
        raise ValueError(f"{path}: expected three integer class codes per marker")
    plt = _pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, name, values in zip(axes, ("J_parallel", "Topology", "Fractal"), classes.T):
        labels, counts = np.unique(values.astype(int), return_counts=True)
        bars = ax.bar(labels.astype(str), counts)
        ax.bar_label(bars, padding=2)
        ax.margins(y=0.15)
        ax.set(title=name, xlabel="Class code (0 = unclassified)")
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Markers")
    fig.suptitle("SIMPLE orbit-class statistics")
    return _finish(fig, output, show)


def plot_wall(value, output=None, color_by="energy", show=None):
    """R-Z wall footprint from authoritative Cartesian hits, with E/E0 or incidence."""
    if color_by not in ("energy", "incidence"):
        raise ValueError("color_by must be 'energy' or 'incidence'")
    diagnostic = "zend" if color_by == "energy" else "wall_hit_cos_incidence"
    path = _results_path(value)
    variables, units = _netcdf(path, ("wall_hit", "wall_hit_cart", diagnostic))
    hit = variables["wall_hit"] == 1
    positions = variables["wall_hit_cart"]
    if hit.ndim != 1 or positions.shape != (len(hit), 3):
        raise ValueError(
            f"{path}: expected one wall flag and Cartesian triple per marker"
        )
    if not np.any(hit):
        raise ValueError(f"{path} contains no wall-hit markers")
    positions = positions[hit]
    color = variables[diagnostic]
    if color_by == "energy":
        if color.shape != (len(hit), 5):
            raise ValueError(f"{path}: expected five zend components per marker")
        color = color[hit, 3] ** 2
        label = r"Impact energy $E/E_0 = (v/v_0)^2$"
    else:
        if color.shape != hit.shape:
            raise ValueError(f"{path}: incidence and wall flag counts differ")
        color = color[hit]
        label = r"$|\cos(\alpha)|$ (0 = grazing, 1 = normal)"
    if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(color)):
        raise ValueError(f"{path} contains missing/nonfinite wall-hit diagnostics")
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(7, 5))
    scatter = ax.scatter(
        np.hypot(positions[:, 0], positions[:, 1]),
        positions[:, 2],
        c=color,
        s=16,
        cmap="viridis",
    )
    fig.colorbar(scatter, ax=ax, label=label)
    ax.set(
        xlabel=f"R / {units}",
        ylabel=f"Z / {units}",
        aspect="equal",
        title=f"SIMPLE wall-hit footprint ({np.count_nonzero(hit)} markers)",
    )
    ax.grid(alpha=0.25)
    return _finish(fig, output, show)


def plot_poincare(directory, output=None, show=None):
    """R-Z sections from spectre_poincare.x; colour is the volume-local s coordinate."""
    directory = Path(directory)
    files = sorted(directory.glob("poincare_vol*.dat"))
    if not files:
        raise FileNotFoundError(f"no poincare_vol<N>.dat files found in {directory}")
    sections = [_table(path, 6) for path in files]
    for path, data in zip(files, sections):
        if not np.all(np.isfinite(data[:, 2:6])):
            raise ValueError(f"{path} contains nonfinite Poincare coordinates")
    plt = _pyplot()
    from matplotlib.colors import Normalize

    norm = Normalize(
        min(d[:, 2].min() for d in sections), max(d[:, 2].max() for d in sections)
    )
    fig, ax = plt.subplots(figsize=(7, 6))
    for path, data in zip(files, sections):
        scatter = ax.scatter(
            data[:, 4],
            data[:, 5],
            c=data[:, 2],
            norm=norm,
            cmap="viridis",
            s=2,
            label=path.stem,
        )
    ax.set(
        xlabel="R (equilibrium units)",
        ylabel="Z (equilibrium units)",
        aspect="equal",
        title="SIMPLE Poincare sections",
    )
    if len(files) <= 12:
        ax.legend(loc="best", fontsize="small")
    fig.colorbar(scatter, ax=ax, label="s (volume-local)")
    return _finish(fig, output, show)


def _parser():
    parser = argparse.ArgumentParser(
        prog="simple-plot",
        description="Plot existing SIMPLE output files; never runs SIMPLE.",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    losses = sub.add_parser(
        "losses",
        help="confinement and loss fraction, with passing/trapped split",
        description=(
            "Pointwise Wilson binomial band on lost fraction among resolved markers. "
            "For p=k/n (k lost, n resolved markers), z=Normal quantile at (1+confidence)/2: "
            "(p+z^2/(2n) +/- z*sqrt(p*(1-p)/n+z^2/(4n^2)))/(1+z^2/n). "
            "Assumes independent, equally weighted markers; excludes numerical failures."
        ),
    )
    losses.add_argument(
        "directory", type=Path, help="directory containing confined_fraction.dat"
    )
    losses.add_argument(
        "--confidence",
        type=float,
        default=0.95,
        metavar="LEVEL",
        help="confidence level between 0 and 1 (default: 0.95)",
    )
    classes = sub.add_parser("classes", help="orbit-class code counts")
    classes.add_argument("path", type=Path, help="output directory or results.nc")
    wall = sub.add_parser("wall", help="R-Z wall-hit footprint from results.nc")
    wall.add_argument("path", type=Path, help="output directory or results.nc")
    wall.add_argument(
        "--color",
        choices=("energy", "incidence"),
        default="energy",
        help="normalized impact energy E/E0 or absolute incidence cosine",
    )
    poincare = sub.add_parser("poincare", help="R-Z sections from poincare_vol<N>.dat")
    poincare.add_argument(
        "directory", type=Path, help="directory containing Poincare files"
    )
    for command in (losses, classes, wall, poincare):
        command.add_argument(
            "-o",
            "--output",
            type=Path,
            help="save a figure (e.g. PNG/PDF); otherwise display interactively",
        )
    return parser


def main(argv=None):
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "losses":
            plot_losses(args.directory, args.output, args.confidence)
        elif args.command == "classes":
            plot_classes(args.path, args.output)
        elif args.command == "wall":
            plot_wall(args.path, args.output, args.color)
        else:
            plot_poincare(args.directory, args.output)
    except (OSError, ImportError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
