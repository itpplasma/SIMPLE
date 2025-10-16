"""Sphinx configuration for the SIMPLE Python API documentation."""

from __future__ import annotations

import importlib.util
import os
import sys
from datetime import datetime
from pathlib import Path

# -- Path setup ----------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python"

sys.path.insert(0, str(PYTHON_SRC))


# -- Project information -------------------------------------------------------

project = "SIMPLE Python API"
author = "SIMPLE Developers"
current_year = datetime.utcnow().year
copyright = f"{current_year}, SIMPLE Developers"
release = os.environ.get("SIMPLE_VERSION", "latest")


# -- General configuration -----------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_mock_imports = ["pysimple"]
autodoc_typehints = "description"
autodoc_type_aliases = {
    "ParticleBatch": "simple.particles.ParticleBatch",
    "BatchResults": "simple.results.BatchResults",
}

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True


# -- HTML output ---------------------------------------------------------------

html_theme = "furo"
html_title = "SIMPLE Python API"
html_static_path = ["_static"]


def skip_private(app, what, name, obj, skip, options):  # noqa: D401
    """Skip private members (starting with `_`) from the documentation."""
    if name.startswith("_"):
        return True
    return skip


def setup(app):  # noqa: D401
    """Register custom Sphinx hooks."""
    app.connect("autodoc-skip-member", skip_private)
