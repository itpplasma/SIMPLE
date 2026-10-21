"""Local small-repository tests for the staged artifact gate."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

GATE = Path(__file__).resolve().parents[2] / "tools/check_staged_artifacts.py"


class StagedArtifactsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        subprocess.run(["git", "init", "-q", str(self.repo)], check=True)

    def add(self, name, contents=b"x"):
        path = self.repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)
        subprocess.run(["git", "add", "--", name], cwd=self.repo, check=True)

    def check(self, override=False):
        env = os.environ.copy()
        env.pop("SIMPLE_ALLOW_BULK_ARTIFACTS", None)
        if override:
            env["SIMPLE_ALLOW_BULK_ARTIFACTS"] = "1"
        return subprocess.run([sys.executable, str(GATE)], cwd=self.repo, env=env,
                              capture_output=True, text=True)

    def test_small_fixture_allowed(self):
        self.add("test/tests/data/fixture.h5", b"small")
        self.assertEqual(self.check().returncode, 0)

    def test_many_small_archives_blocked_unless_explicit_override(self):
        for i in range(25):
            self.add(f"concept/results/row{i:03}.csv.gz")
        self.assertEqual(self.check().returncode, 1)
        self.assertEqual(self.check(override=True).returncode, 0)

    def test_large_aggregate_blocked(self):
        self.add("concept/results/trace.npz", b"x" * (4 * 1024 * 1024 + 1))
        self.assertEqual(self.check().returncode, 1)


if __name__ == "__main__":
    unittest.main()
