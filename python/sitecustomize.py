from __future__ import annotations

import sys
from pathlib import Path


def _prepend_sys_path(path: Path) -> None:
    path_str = str(path)
    try:
        sys.path.remove(path_str)
    except ValueError:
        pass
    sys.path.insert(0, path_str)


_REPO_ROOT = Path(__file__).resolve().parents[1]

_prepend_sys_path(_REPO_ROOT / "build" / "python")
_prepend_sys_path(_REPO_ROOT / "build")
_prepend_sys_path(_REPO_ROOT / "python")
