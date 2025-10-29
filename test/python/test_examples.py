#!/usr/bin/env python3
"""
Test that all example scripts work correctly.

This ensures documentation examples stay valid as single source of truth.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

pytest.importorskip("pysimple", reason="pysimple module not available")

examples_dir = Path(__file__).resolve().parents[2] / "examples"
sys.path.insert(0, str(examples_dir))

import simple_api
import classify_fast
import orbit_macro


class TestExamples:
    """Test that all documented examples execute without errors."""

    def test_simple_api_trace_example(self, vmec_file: str):
        """Test simple_api.py trace example."""
        simple_api.trace_example(vmec_file, n_particles=4)

    def test_simple_api_classify_example(self, vmec_file: str):
        """Test simple_api.py classify example."""
        simple_api.classify_example(vmec_file, n_particles=4)

    def test_classify_fast_example(self, vmec_file: str):
        """Test classify_fast.py example."""
        result = classify_fast.classify_fast_example(vmec_file, n_particles=4)
        assert 'passing' in result
        assert 'jpar' in result

    def test_orbit_macro_example(self, vmec_file: str):
        """Test orbit_macro.py example."""
        result = orbit_macro.trace_macrostep_example(vmec_file)
        assert 'trajectory' in result
        assert 'times' in result
        assert result['trajectory'].shape[0] == 5
