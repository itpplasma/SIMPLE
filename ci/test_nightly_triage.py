"""Guards for the nightly triage logic, with no GitHub or git access."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nightly_triage as nt  # noqa: E402


class TriageTests(unittest.TestCase):
    def test_config_names_every_job_with_a_command(self):
        jobs = nt.load_jobs()
        self.assertTrue(jobs)
        workflow = (Path(__file__).resolve().parents[1] / ".github" / "workflows"
                    / nt.WORKFLOW).read_text()
        for name in jobs:
            self.assertIn(f'nightly_triage.py run "{name}"', workflow, name)

    def test_failing_jobs(self):
        configured = {"unit": "true", "slow (Release)": "true"}
        jobs = [{"name": "unit", "conclusion": "failure"},
                {"name": "nightly / slow (Release)", "conclusion": "timed_out"},
                {"name": "Bisect failing jobs and report", "conclusion": "failure"},
                {"name": "other", "conclusion": "success"}]
        self.assertEqual(nt.failing_jobs(jobs, configured), ["slow (Release)", "unit"])
        jobs[0]["conclusion"] = "success"
        self.assertEqual(nt.failing_jobs(jobs, configured), ["slow (Release)"])

    def test_last_green_skips_this_run_and_red_runs(self):
        runs = [{"databaseId": 9, "conclusion": None, "headSha": "now"},
                {"databaseId": 8, "conclusion": "failure", "headSha": "red"},
                {"databaseId": 7, "conclusion": "success", "headSha": "green"}]
        self.assertEqual(nt.last_green(runs, 9)["headSha"], "green")
        self.assertIsNone(nt.last_green(runs[:2], 9))

    def test_report_records_jobs_without_mentions(self):
        body = nt.failure_report(
            branch="main", bad="b" * 40, good="a" * 40, commits=3,
            culprit={"sha": "c" * 40, "subject": "Break the slow tier"}, pull=42,
            failed=["slow (Release)", "unit"], run_url="https://x/run",
            job_urls={"unit": "https://x/j"}, server="https://github.com", repo="o/r")
        self.assertIn("slow (Release), [unit](https://x/j)", body)
        self.assertIn("First bad commit: " + "c" * 40 + " Break the slow tier from #42", body)
        self.assertIn("compare/" + "a" * 40 + "..." + "b" * 40, body)
        self.assertNotIn("@", body)
        self.assertEqual(nt.issue_jobs([body, "a comment"]), {"slow (Release)", "unit"})

    def test_step_finds_a_bad_commit(self):
        with tempfile.TemporaryDirectory() as tmp:
            def git(*args):
                subprocess.run(["git", "-C", tmp, *args], check=True, capture_output=True)
            git("init", "-q")
            Path(tmp, "ok").write_text("1")
            git("add", "ok")
            git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "good")
            step = [sys.executable, nt.__file__, "step", "--command", "grep -q 1 ok"]
            self.assertEqual(subprocess.run(step, cwd=tmp).returncode, 0)
            Path(tmp, "ok").write_text("0")
            self.assertEqual(subprocess.run(step, cwd=tmp).returncode, 1)


if __name__ == "__main__":
    unittest.main()
