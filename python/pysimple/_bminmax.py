"""bminmax cache lifecycle helpers."""

from __future__ import annotations


def needs_bminmax_cache(num_surf: int, ntcut: int, class_plot: bool) -> bool:
    """Match the Fortran cache-reader predicate."""
    if int(ntcut) > 0 or bool(class_plot):
        return int(num_surf) > 1
    return int(num_surf) != 1
