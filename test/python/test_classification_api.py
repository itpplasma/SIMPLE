#!/usr/bin/env python3
"""Tests for the fast classification Python API."""

from __future__ import annotations

import os
import math
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")

import pysimple

REPO_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_EXE = REPO_ROOT / "build" / "simple.x"


@pytest.mark.usefixtures("vmec_file")
def test_classify_parallel_api(vmec_file: str) -> None:
    # Explicitly set ntestpart to avoid interference
    pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, tcut=0.1, ntestpart=6)
    particles = pysimple.sample_surface(6, s=0.35)

    results = pysimple.classify_parallel(particles)

    assert results['jpar'].shape == (6,)
    assert results['topology'].shape == (6,)
    assert results['fractal'].shape == (6,)
    assert results['passing'].shape == (6,)


def test_classification_values(vmec_file: str) -> None:
    # Explicitly set ntestpart to avoid interference
    pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, tcut=0.1, ntestpart=8)
    particles = pysimple.sample_surface(8, s=0.4)

    results = pysimple.classify_parallel(particles)

    assert results['jpar'].min() >= 0
    assert results['jpar'].max() <= 2
    assert results['topology'].min() >= 0
    assert results['topology'].max() <= 2
    assert results['fractal'].min() >= 0
    assert results['fractal'].max() <= 3
