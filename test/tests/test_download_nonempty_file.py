#!/usr/bin/env python3
"""Exercise CMake downloads with deterministic local files and failed transfers."""
from pathlib import Path
import subprocess
import sys
import tempfile

helper = Path(sys.argv[1]).resolve()
cmake = sys.argv[2]
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    source = root / "source.dat"
    source.write_bytes(b"synthetic test fixture\n")
    destination = root / "download.dat"
    script = root / "download.cmake"

    def download(url):
        script.write_text(
            f'include("{helper}")\n'
            f'download_nonempty_file("{url}" "{destination}")\n'
        )
        return subprocess.run([cmake, "-P", script], text=True,
                              capture_output=True, timeout=10)

    destination.touch()
    completed = download(source.as_uri())
    assert completed.returncode == 0, completed.stderr
    assert destination.read_bytes() == source.read_bytes(), "empty cache was reused"
    source.write_bytes(b"changed source\n")
    completed = download(source.as_uri())
    assert completed.returncode == 0, completed.stderr
    assert destination.read_bytes() == b"synthetic test fixture\n", "cache overwritten"
    destination.unlink()
    completed = download((root / "absent.dat").as_uri())
    assert completed.returncode != 0, "failed transfer was silently accepted"
    assert not destination.exists(), "failed transfer left a poisoned cache"
    source.write_bytes(b"")
    completed = download(source.as_uri())
    assert completed.returncode != 0, "empty transfer was silently accepted"
    assert not destination.exists(), "empty transfer left a poisoned cache"
print("CMake download cache and failure handling passed")
