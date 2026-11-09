"""
Pytest configuration and fixtures for SIMPLE tests.

Provides shared test data management, including automatic download
of VMEC equilibrium file.
"""

import sys
from pathlib import Path
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))
sys.path.insert(0, str(_REPO_ROOT / "build"))
sys.path.insert(0, str(_REPO_ROOT / "build" / "python"))

from test_data_fetch import WOUT_QA_SHA256, WOUT_QA_URL, fetch

# Test data location - single source of truth for all tests
TEST_DATA_DIR = Path(__file__).parent / "test_data"
WOUT_FILE = TEST_DATA_DIR / "wout.nc"


def download_wout_file():
    """Download VMEC equilibrium file if not present or altered."""
    fetch(WOUT_QA_URL, WOUT_FILE, WOUT_QA_SHA256)


@pytest.fixture(scope="session")
def vmec_file():
    """Provide path to VMEC equilibrium file, downloading if necessary.

    This fixture is session-scoped, so the file is only downloaded once
    per test session and shared across all tests.

    Returns:
        Path: Absolute path to wout.nc file
    """
    download_wout_file()
    return str(WOUT_FILE.absolute())


@pytest.fixture(scope="session", autouse=True)
def setup_test_data():
    """Auto-download test data at the start of test session."""
    download_wout_file()
