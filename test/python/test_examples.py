#!/usr/bin/env python3
"""
Test that all example scripts execute without errors.

This ensures documentation examples stay valid as single source of truth.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")


class TestExamples:
    """Test that all documented examples execute without errors."""

    def test_simple_api_example(self, vmec_file: str):
        """Test simple_api.py example executes."""
        examples_dir = Path(__file__).resolve().parents[2] / "examples"
        result = subprocess.run(
            [sys.executable, str(examples_dir / "simple_api.py")],
            capture_output=True,
            text=True,
            cwd=str(examples_dir),
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    def test_classify_fast_example(self, vmec_file: str):
        """Test classify_fast.py example executes."""
        examples_dir = Path(__file__).resolve().parents[2] / "examples"
        result = subprocess.run(
            [sys.executable, str(examples_dir / "classify_fast.py")],
            capture_output=True,
            text=True,
            cwd=str(examples_dir),
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    def test_classify_fractal_example(self, vmec_file: str):
        """Test classify_fractal.py example executes."""
        examples_dir = Path(__file__).resolve().parents[2] / "examples"
        result = subprocess.run(
            [sys.executable, str(examples_dir / "classify_fractal.py")],
            capture_output=True,
            text=True,
            cwd=str(examples_dir),
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"

    def test_orbit_macro_example(self, vmec_file: str):
        """Test orbit_macro.py example executes."""
        examples_dir = Path(__file__).resolve().parents[2] / "examples"
        result = subprocess.run(
            [sys.executable, str(examples_dir / "orbit_macro.py")],
            capture_output=True,
            text=True,
            cwd=str(examples_dir),
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"
