"""bminmax cache lifecycle helpers."""

from __future__ import annotations


def classification_enabled(ntcut: int, class_plot: bool, fast_class: bool) -> bool:
    """Match the Fortran params.classification_enabled predicate."""
    return int(ntcut) > 0 or bool(class_plot) or bool(fast_class)


def needs_bminmax_cache(
    num_surf: int, ntcut: int, class_plot: bool, fast_class: bool = False
) -> bool:
    """Match the Fortran cache-reader predicate."""
    if classification_enabled(ntcut, class_plot, fast_class):
        return int(num_surf) > 1
    return int(num_surf) != 1
