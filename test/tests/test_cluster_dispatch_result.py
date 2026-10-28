"""Run the cluster workflow step against dispatcher completion outcomes."""

from pathlib import Path
import os
import subprocess
import tempfile
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[2]


def workflow_step():
    lines = (ROOT / ".github/workflows/main.yml").read_text().splitlines()
    start = next(
        i for i, line in enumerate(lines)
        if line.strip() == "- name: Build and test on sCluster"
    )
    start = next(i for i in range(start, len(lines)) if lines[i].strip() == "run: |")
    indent = len(lines[start]) - len(lines[start].lstrip())
    body = []
    for line in lines[start + 1:]:
        if line.strip() and len(line) - len(line.lstrip()) <= indent:
            break
        body.append(line)
    return textwrap.dedent("\n".join(body))


class ClusterDispatchResult(unittest.TestCase):
    def run_step(self, outcome, event="push", full="false", stale=False):
        with tempfile.TemporaryDirectory(prefix="cluster-result-") as tmp:
            root = Path(tmp)
            for directory in ("ci", "bin", "sci-out"):
                (root / directory).mkdir()
            bundle = root / "ci/bundle.sh"
            bundle.write_text("#!/bin/bash\nexit 0\n")
            bundle.chmod(0o755)
            dispatch = root / "bin/sci-dispatch"
            dispatch.write_text(
                '#!/bin/bash\nset -eu\n'
                'case $DISPATCH_RESULT in\n'
                'pass) printf "Correctness gate: exit 0 (authoritative)\\n" '
                '> "$SCI_OUT/full-status.txt" ;;\n'
                'fail) printf "Correctness gate: exit 1 (authoritative)\\n'
                'Coverage worker: exit 0 (reporting only)\\n" '
                '> "$SCI_OUT/full-status.txt" ;;\n'
                'missing) ;;\nesac\nexit 0\n'
            )
            dispatch.chmod(0o755)
            if stale:
                (root / "sci-out/full-status.txt").write_text(
                    "Correctness gate: exit 0 (authoritative)\n"
                )
            env = dict(os.environ, EVENT=event, FULL_INPUT=full,
                       RUNNER_TEMP=tmp, SCI_OUT=str(root / "sci-out"),
                       DISPATCH_RESULT=outcome,
                       PATH=str(root / "bin") + os.pathsep + os.environ["PATH"])
            return subprocess.run(
                ["bash", "-e", "-c", workflow_step()], cwd=root, env=env,
                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=10,
            )

    def test_full_requires_completed_results(self):
        for stale in (False, True):
            with self.subTest(stale=stale):
                result = self.run_step("missing", stale=stale)
                self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_reporting_success_cannot_override_failed_gate(self):
        result = self.run_step("fail")
        self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_completed_full_gate_passes(self):
        result = self.run_step("pass")
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_manual_smoke_does_not_require_full_results(self):
        result = self.run_step("missing", event="workflow_dispatch")
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_merge_queue_runs_the_fast_tier(self):
        result = self.run_step("missing", event="merge_group")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("tier: fast", result.stdout)

    def test_manual_full_requires_completed_results(self):
        result = self.run_step("missing", event="workflow_dispatch", full="true")
        self.assertNotEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
