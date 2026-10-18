"""Shared matplotlib style for SIMPLE documentation figures.

Import this module from figure scripts in this folder and call ``use()``
before plotting. Categorical colours follow the Okabe-Ito palette, which is
distinguishable under the common colour-vision deficiencies; always pair a
colour with a second cue (line style, marker or direct label). Continuous data
use perceptually uniform maps (``SEQUENTIAL`` or ``SEQUENTIAL_ALT``).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

# Okabe & Ito (2008), "Color Universal Design".
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}

# Semantic roles used across SIMPLE figures.
COLORS = {
    "total": OKABE_ITO["black"],
    "passing": OKABE_ITO["blue"],
    "trapped": OKABE_ITO["vermillion"],
    "reference": OKABE_ITO["green"],
    "accent": OKABE_ITO["orange"],
    "neutral": "#6E6E6E",
    "surface": "#B9C2CC",
}

LINESTYLES = {"total": "-", "passing": "--", "trapped": ":"}

SEQUENTIAL = "viridis"
SEQUENTIAL_ALT = "cividis"

# Text and grid colours chosen to stay legible when the figure is placed on an
# opaque white card in either a light or a dark page theme.
TEXT = "#222222"
GRID = "#D9D9D9"

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "axes.edgecolor": TEXT,
    "axes.labelcolor": TEXT,
    "text.color": TEXT,
    "xtick.color": TEXT,
    "ytick.color": TEXT,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.prop_cycle": mpl.cycler(
        color=[
            OKABE_ITO["blue"],
            OKABE_ITO["vermillion"],
            OKABE_ITO["green"],
            OKABE_ITO["orange"],
            OKABE_ITO["purple"],
            OKABE_ITO["sky"],
            OKABE_ITO["black"],
        ]
    ),
    "lines.linewidth": 1.6,
    "image.cmap": SEQUENTIAL,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    "svg.fonttype": "none",  # keep text as text in SVG output
    "pdf.fonttype": 42,
}


def use() -> None:
    """Apply the SIMPLE figure style to the current matplotlib session."""
    mpl.rcParams.update(RC)


def save(fig, path: str | Path, **kwargs) -> Path:
    """Save ``fig`` with the shared defaults and return the output path.

    PNG output gets an opaque white background so that dark page themes do not
    swallow dark text; pass ``transparent=True`` to override.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".png":
        kwargs.setdefault("metadata", {"Software": None})
        if not kwargs.get("transparent", False):
            kwargs.setdefault("facecolor", "white")
    fig.savefig(path, **kwargs)
    return path
