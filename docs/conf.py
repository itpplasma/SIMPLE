"""Sphinx configuration for the SIMPLE documentation."""

from __future__ import annotations

import importlib.util
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# -- Path setup ----------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python"
PYTHON_SRC_STR = str(PYTHON_SRC)
if PYTHON_SRC_STR not in sys.path:
    sys.path.insert(0, PYTHON_SRC_STR)


# -- Project information -------------------------------------------------------

project = "SIMPLE"
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
    "sphinx.ext.mathjax",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 6

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build", "Thumbs.db", ".DS_Store"]

DOCS_DIR = Path(__file__).resolve().parent
FIGURE_DIR = DOCS_DIR / "_static" / "figures"
PLACEHOLDER_MARKER = ".. placeholder-page"

# A generated input reference may arrive as input-reference.md. Then the
# committed placeholder input-reference.rst must not shadow it.
_generated_md = DOCS_DIR / "input-reference.md"
_placeholder_rst = DOCS_DIR / "input-reference.rst"
if _generated_md.exists() and _placeholder_rst.exists():
    if PLACEHOLDER_MARKER in _placeholder_rst.read_text(encoding="utf-8"):
        exclude_patterns.append("input-reference.rst")

autodoc_class_signature = "separated"
autodoc_member_order = "bysource"
autodoc_mock_imports = ["simple_backend", "_simple_main", "_simple_backend"]
autodoc_typehints = "description"

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True


# -- HTML output ---------------------------------------------------------------

html_theme = "furo"
html_title = "SIMPLE documentation"
html_static_path = ["_static"]


def skip_private(app, what, name, obj, skip, options):  # noqa: D401
    """Skip private members (starting with `_`) from the documentation."""
    if name.startswith("_"):
        return True
    return skip


_FIGURE_RE = re.compile(
    r"^(?P<indent>[ \t]*)\.\. figure:: /_static/figures/(?P<name>\S+)\n"
    r"(?P<body>(?:(?P=indent)[ \t]+\S.*\n|[ \t]*\n)*)",
    re.MULTILINE,
)


def placeholder_missing_figures(app, docname, source):  # noqa: D401
    """Replace figures under _static/figures that are not yet committed.

    Figures such as the pipeline diagram are produced separately. Until a
    file exists, its figure directive becomes a note, so that the build with
    warnings as errors does not fail.
    """

    def replace(match: re.Match[str]) -> str:
        if (FIGURE_DIR / match["name"]).exists():
            return match[0]
        indent = match["indent"]
        return (
            f"{indent}.. note::\n\n"
            f"{indent}   Figure ``{match['name']}`` is not yet available.\n\n"
        )

    source[0] = _FIGURE_RE.sub(replace, source[0])


def setup(app):  # noqa: D401
    """Register custom Sphinx hooks."""
    app.connect("autodoc-skip-member", skip_private)
    app.connect("source-read", placeholder_missing_figures)
