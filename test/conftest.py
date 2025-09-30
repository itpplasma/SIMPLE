"""
Pytest configuration and fixtures for SIMPLE tests.

Provides shared test data management, including automatic download
of VMEC equilibrium file.
"""

import pytest
from pathlib import Path
import urllib.request
import hashlib

# Test data location - single source of truth for all tests
TEST_DATA_DIR = Path(__file__).parent / "test_data"
WOUT_FILE = TEST_DATA_DIR / "wout.nc"
WOUT_URL = "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

# Expected SHA256 hash for verification (optional but recommended)
# To generate: sha256sum wout.nc
# WOUT_SHA256 = "..."  # TODO: Add hash for verification


def download_wout_file():
    """Download VMEC equilibrium file if not present."""
    if WOUT_FILE.exists():
        return

    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading test VMEC file to {WOUT_FILE}...")
    try:
        urllib.request.urlretrieve(WOUT_URL, WOUT_FILE)
        print(f"Downloaded {WOUT_FILE.stat().st_size} bytes")
    except Exception as e:
        raise RuntimeError(f"Failed to download VMEC test file: {e}")


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