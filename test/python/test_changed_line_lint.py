"""Exercise the changed-line gate with real Git diffs, Ruff and fprettify."""

import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "tools/check_changed_lines.py"
SPEC = importlib.util.spec_from_file_location("changed_line_lint", CHECKER)
lint = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lint)


class SourceCoordinates(unittest.TestCase):
    def test_formatter_wrapping_does_not_shift_the_selected_source_line(self):
        source = 'legacy = ["' + "x" * 100 + '"]\n\nchanged=2\n'
        formatted = 'legacy = [\n    "' + "x" * 100 + '",\n]\n\nchanged = 2\n'
        self.assertIn(3, lint.formatting_lines(source, formatted))
        self.assertNotIn(2, lint.formatting_lines(source, formatted))

    def test_formatter_deletions_use_source_coordinates(self):
        self.assertEqual(
            lint.formatting_lines("x = 1\n\n\n\ny = 2\n", "x = 1\n\n\ny = 2\n"), {4}
        )


@unittest.skipUnless(
    shutil.which("ruff") and shutil.which("fprettify"),
    "install the lint tools to run integration checks",
)
class ChangedLineGate(unittest.TestCase):
    def setUp(self):
        (ROOT / "build").mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=ROOT / "build")
        self.addCleanup(self.temporary.cleanup)
        self.repo = Path(self.temporary.name)
        self.git("init", "-q")
        self.git("config", "user.name", "Lint test")
        self.git("config", "user.email", "lint@example.invalid")
        self.git("config", "commit.gpgsign", "false")
        shutil.copy(ROOT / ".fprettify", self.repo)
        (self.repo / "tools").mkdir()
        self.python = self.repo / "tools" / "legacy ü space.py"
        self.fortran = self.repo / "legacy.F90"
        self.python.write_text(
            "legacy=1\nold_error = unknown_legacy_name\n\nvalue = 2\n"
        )
        self.fortran.write_text(
            "program example\ninteger::legacy\n    legacy=1\nend program example\n"
        )
        self.git("add", ".")
        self.git("commit", "-qm", "legacy baseline")
        self.base = self.git("rev-parse", "HEAD").stdout.strip()

    def git(self, *args):
        return subprocess.run(
            ["git", *args], cwd=self.repo, text=True, capture_output=True, check=True
        )

    def check(self, expected):
        before = (self.python.read_bytes(), self.fortran.read_bytes())
        result = subprocess.run(
            [sys.executable, str(CHECKER), "--base", self.base],
            cwd=self.repo,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        self.assertEqual(before, (self.python.read_bytes(), self.fortran.read_bytes()))
        return result.stdout

    def test_good_changed_lines_ignore_legacy_format_and_diagnostics(self):
        self.python.write_text(
            self.python.read_text().replace("value = 2", "value = 3")
        )
        self.fortran.write_text(
            self.fortran.read_text().replace("legacy=1", "legacy = 2")
        )
        self.check(0)

    def test_changed_python_format_after_legacy_line_wrap(self):
        baseline = 'legacy = ["' + "x" * 100 + '"]\n\nvalue = 2\n'
        self.python.write_text(baseline)
        self.git("add", ".")
        self.git("commit", "-qm", "long legacy line")
        self.base = self.git("rev-parse", "HEAD").stdout.strip()
        self.python.write_text(baseline.replace("value = 2", "value = 3"))
        self.check(0)
        self.python.write_text(baseline.replace("value = 2", "value=3"))
        self.git("add", ".")
        self.git("commit", "-qm", "bad formatting")
        self.assertIn("legacy ü space.py:3: formatting differs", self.check(1))

    def test_changed_ruff_diagnostic_is_rejected(self):
        self.python.write_text(
            self.python.read_text().replace("value = 2", "value = missing")
        )
        self.assertIn("F821", self.check(1))

    def test_bad_changed_fortran_format_is_rejected(self):
        self.fortran.write_text(
            self.fortran.read_text().replace("legacy=1", "legacy=2")
        )
        self.git("add", ".")  # Also exercise staged changes.
        self.assertIn("legacy.F90:3: formatting differs", self.check(1))

    def test_deletion_of_legacy_line_does_not_check_its_neighbor(self):
        self.python.write_text(self.python.read_text().replace("value = 2\n", ""))
        self.check(0)

    def test_rename_keeps_untouched_legacy_lines_out_of_the_gate(self):
        destination = self.repo / "tools" / "renamed.py"
        self.git(
            "mv",
            str(self.python.relative_to(self.repo)),
            str(destination.relative_to(self.repo)),
        )
        self.python = destination
        self.check(0)
        self.python.write_text(
            self.python.read_text().replace("value = 2", "value = 3")
        )
        self.check(0)
        self.python.write_text(self.python.read_text().replace("value = 3", "value=3"))
        self.assertIn("renamed.py:4: formatting differs", self.check(1))

    def test_new_file_and_invalid_base(self):
        (self.repo / "tools" / "new.py").write_text("value=1\n")
        self.git("add", ".")
        self.assertIn("new.py:1", self.check(1))
        self.base = "nonexistent-revision"
        self.check(2)


if __name__ == "__main__":
    unittest.main()
