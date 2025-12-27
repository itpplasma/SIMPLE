"""
Pytest configuration and fixtures for SIMPLE tests.

Provides shared test data management, including automatic download
of VMEC equilibrium file.
"""

import sys
from pathlib import Path
import urllib.request
import hashlib
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))
sys.path.insert(0, str(_REPO_ROOT / "build"))
sys.path.insert(0, str(_REPO_ROOT / "build" / "python"))

# Test data location - single source of truth for all tests
TEST_DATA_DIR = Path(__file__).parent / "test_data"
WOUT_FILE = TEST_DATA_DIR / "wout.nc"
WOUT_URL = "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

WOUT_SHA256 = "333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771"


def _sha256sum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ensure_wout_is_valid() -> None:
    if not WOUT_FILE.exists():
        return
    if _sha256sum(WOUT_FILE) == WOUT_SHA256:
        return
    WOUT_FILE.unlink()


def download_wout_file():
    """Download VMEC equilibrium file if not present."""
    _ensure_wout_is_valid()
    if WOUT_FILE.exists():
        return

    TEST_DATA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Downloading test VMEC file to {WOUT_FILE}...")
    try:
        urllib.request.urlretrieve(WOUT_URL, WOUT_FILE)
        if _sha256sum(WOUT_FILE) != WOUT_SHA256:
            WOUT_FILE.unlink()
            raise RuntimeError("Downloaded VMEC test file has unexpected SHA256 hash")
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
