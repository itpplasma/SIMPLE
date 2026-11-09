"""Fetch pinned test data: retried on transient errors, checked by SHA-256.

URLs and hashes match test/tests/CMakeLists.txt and ci/bundle.sh.
"""
import hashlib
import time
import urllib.error
import urllib.request
from pathlib import Path

SIMSOPT_TEST_FILES = (
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/"
    "912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files")
WOUT_QA_URL = (SIMSOPT_TEST_FILES
               + "/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc")
WOUT_QA_SHA256 = "333ab8cef64c4b4f406e76209266fc3bf4d7976e9c28c162983c2120e112e771"
# SHA-256 of every pinned file, by its name at the source.
PINNED_SHA256 = {
    "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc": WOUT_QA_SHA256,
    "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc":
        "249c397d67357949b2ab4690145224736d8388f5332f23c67c9cb0bf2aa9931e",
    "wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc":
        "36926d748099bc7abc138d82bdd4b447d2579f08f50189f27d441cf96c1b51e0",
    "coils.c09r00": "3c429da06f4c062887a497a16e2d2bd10f0ecb0b8858c252698631f3853da428",
    "wout_ncsx_c09r00_fixed.nc":
        "1005d81729366118ed51926c51b99b08bea21cfa77bb6f2874774632fb0a1e99",
}


def sha256sum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(url, path, sha256, attempts=5):
    """Make path hold the file at url with the given SHA-256."""
    path = Path(path)
    if path.exists() and sha256sum(path) == sha256:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                part.write_bytes(response.read())
            actual = sha256sum(part)
            if actual == sha256:
                part.replace(path)
                return path
            error = f"SHA-256 {actual}, expected {sha256}"
        except (urllib.error.URLError, OSError) as exc:
            error = str(exc)
        if attempt < attempts:
            time.sleep(10 * attempt)
    part.unlink(missing_ok=True)
    raise RuntimeError(f"cannot fetch {url}: {error}")


def fetch_pinned(url, path):
    """fetch with the SHA-256 pinned for the file name in url."""
    return fetch(url, path, PINNED_SHA256[url.rsplit("/", 1)[-1]])
