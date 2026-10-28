#!/usr/bin/env python3
"""Triage a failed nightly run: bisect the failing jobs and file an issue.

The jobs and the command each one runs live in nightly-jobs.json next to this
script, so the nightly workflow and the bisection run the same command.

Subcommands, run inside GitHub Actions (GH_TOKEN needs issues: write):

  run NAME  run the command of nightly job NAME
  failure   find the failing jobs and the last green nightly, bisect, report
  step      one `git bisect run` step: run the given commands in a clean tree
  success   close the open nightly-failure issue
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "nightly-jobs.json"
WORKFLOW = "nightly.yml"
LABEL = "nightly-failure"
MARKER = re.compile(r"<!-- nightly-failing-jobs: (.*?) -->")
FAILED = ("failure", "timed_out")


def load_jobs(path=CONFIG):
    with open(path) as handle:
        jobs = json.load(handle)
    if not jobs or not all(isinstance(v, str) and v for v in jobs.values()):
        raise ValueError(f"{path}: expected a non-empty map of job name to command")
    return jobs


def job_key(name, configured):
    """Configured job a GitHub job name belongs to, or None.

    A job called from a reusable workflow is reported as "caller / job".
    """
    if name in configured:
        return name
    tail = name.rsplit(" / ", 1)[-1]
    return tail if tail in configured else None


def failing_jobs(jobs, configured):
    names = set()
    for job in jobs:
        key = job_key(job["name"], configured)
        if key and job.get("conclusion") in FAILED:
            names.add(key)
    return sorted(names)


def last_green(runs, current_run):
    """Newest successful nightly run other than this one."""
    for run in runs:
        if run["databaseId"] != current_run and run.get("conclusion") == "success":
            return run
    return None


def issue_title(branch):
    return f"Nightly failure on {branch}"


def issue_jobs(texts):
    """All failing jobs recorded in an issue body and its comments."""
    names = set()
    for text in texts:
        for match in MARKER.finditer(text or ""):
            names.update(json.loads(match.group(1)))
    return names


def failure_report(*, branch, bad, good, commits, culprit, pull, failed,
                   run_url, job_urls, server, repo, note=None):
    lines = [f"The nightly run failed on `{branch}` at {bad}.", ""]
    lines.append("Failing jobs: " + ", ".join(
        f"[{name}]({job_urls[name]})" if name in job_urls else name
        for name in failed))
    if good:
        lines.append(f"Commit range: {good}..{bad} ({commits} commits since the "
                     f"last green nightly, [compare]({server}/{repo}/compare/{good}...{bad}))")
    else:
        lines.append("Commit range: unknown, there is no earlier green nightly run")
    if culprit:
        via = f" from #{pull}" if pull else ""
        lines.append(f"First bad commit: {culprit['sha']} {culprit['subject']}{via}")
    if note:
        lines.append(note)
    lines.append(f"Run and bisection log: {run_url}")
    lines.append("")
    lines.append(f"<!-- nightly-failing-jobs: {json.dumps(failed)} -->")
    return "\n".join(lines)


# ---- GitHub and git plumbing -------------------------------------------------

def api(path, method="GET", data=None):
    """One GitHub REST call with the job token; returns the decoded JSON."""
    url = path if path.startswith("http") else f"{env('GITHUB_API_URL')}/{path}"
    request = urllib.request.Request(
        url, method=method, data=json.dumps(data).encode() if data is not None else None,
        headers={"Authorization": f"Bearer {env('GH_TOKEN')}",
                 "Accept": "application/vnd.github+json",
                 "X-GitHub-Api-Version": "2022-11-28"})
    with urllib.request.urlopen(request, timeout=60) as response:
        body = response.read()
    return json.loads(body) if body else None


def git(*args, check=True):
    return subprocess.run(["git", *args], check=check, capture_output=True,
                          text=True).stdout.strip()


def env(name):
    return os.environ[name]


def open_issue(repo, branch):
    issues = api(f"repos/{repo}/issues?labels={LABEL}&state=open&per_page=100")
    for issue in issues:
        if issue["title"] == issue_title(branch) and "pull_request" not in issue:
            return issue["number"]
    return None


def cmd_run(args):
    command = load_jobs()[args.name]
    print(f"+ {command}", flush=True)
    return subprocess.run(["bash", "-euo", "pipefail", "-c", command]).returncode


def cmd_step(args):
    """Exit 0 for good, 1 for bad; every command must pass."""
    for command in args.command:
        git("clean", "-ffdxq")
        print(f"+ {command}", flush=True)
        if subprocess.run(["bash", "-euo", "pipefail", "-c", command]).returncode:
            return 1
    return 0


def bisect(good, bad, commands):
    """First bad commit in good..bad for the given commands.

    Returns None when the commands also fail at the green commit, since then
    no commit in the range introduced the failure.
    """
    script = Path(os.environ.get("RUNNER_TEMP", "/tmp")) / "nightly_triage.py"
    shutil.copy(__file__, script)  # bisect checks out older versions of this file
    step = [sys.executable, str(script), "step"]
    for command in commands:
        step += ["--command", command]
    git("checkout", "--quiet", good)
    try:
        good_still_passes = subprocess.run(step).returncode == 0
    finally:
        git("clean", "-ffdxq")
        git("checkout", "--quiet", bad)
    if not good_still_passes:
        return None
    git("bisect", "start", bad, good)
    try:
        subprocess.run(["git", "bisect", "run", *step], check=True)
        sha = git("rev-parse", "refs/bisect/bad")
    finally:
        print(git("bisect", "log", check=False))
        git("bisect", "reset", check=False)
    return {"sha": sha, "subject": git("log", "-1", "--format=%s", sha)}


def cmd_failure(args):
    configured = load_jobs()
    repo, run_id = env("GITHUB_REPOSITORY"), int(env("GITHUB_RUN_ID"))
    server, branch, bad = env("GITHUB_SERVER_URL"), env("GITHUB_REF_NAME"), env("GITHUB_SHA")
    jobs = api(f"repos/{repo}/actions/runs/{run_id}/jobs?per_page=100")["jobs"]
    failed = failing_jobs(jobs, configured)
    if args.jobs:
        # Jobs that need another runner are triaged by another triage job.
        failed = [name for name in failed if name in args.jobs.split(",")]
        if not failed:
            print(f"none of {args.jobs} failed")
            return 0
    if not failed:
        raise SystemExit("the nightly failed but no configured job did")
    job_urls = {}
    for job in jobs:
        key = job_key(job["name"], configured)
        if key:
            job_urls[key] = job["html_url"]
    runs = api(f"repos/{repo}/actions/workflows/{WORKFLOW}/runs?per_page=100&"
               + urllib.parse.urlencode({"branch": branch}))["workflow_runs"]
    runs = [{"databaseId": r["id"], "headSha": r["head_sha"],
             "conclusion": r["conclusion"]} for r in runs]
    green = last_green(runs, run_id)
    good = green["headSha"] if green else None
    culprit, pull, note, commits = None, None, None, 0
    if good:
        commits = int(git("rev-list", "--count", f"{good}..{bad}"))
        if good == bad:
            note = "No commit since the last green nightly: the failure is not from a code change."
        elif subprocess.run(["git", "merge-base", "--is-ancestor", good, bad]).returncode:
            note = "The last green commit is not an ancestor; bisection skipped."
        else:
            culprit = bisect(good, bad, [configured[name] for name in failed])
            if culprit:
                pulls = api(f"repos/{repo}/commits/{culprit['sha']}/pulls")
                pull = pulls[0]["number"] if pulls else None
            else:
                note = ("The failing jobs also fail at the last green commit when "
                        "rerun, so the failure is not from a commit in the range.")
    body = failure_report(
        branch=branch, bad=bad, good=good, commits=commits, culprit=culprit,
        pull=pull, failed=failed, run_url=f"{server}/{repo}/actions/runs/{run_id}",
        job_urls=job_urls, server=server, repo=repo, note=note)
    print(body)
    number = open_issue(repo, branch)
    if number:
        api(f"repos/{repo}/issues/{number}/comments", "POST", {"body": body})
    else:
        try:
            api(f"repos/{repo}/labels", "POST", {"name": LABEL, "color": "B60205",
                                                 "description": "The nightly run failed"})
        except urllib.error.HTTPError as error:
            if error.code != 422:  # the label exists
                raise
        api(f"repos/{repo}/issues", "POST",
            {"title": issue_title(branch), "labels": [LABEL], "body": body})


def cmd_success(args):
    repo, branch = env("GITHUB_REPOSITORY"), env("GITHUB_REF_NAME")
    number = open_issue(repo, branch)
    if not number:
        return
    run_url = f"{env('GITHUB_SERVER_URL')}/{repo}/actions/runs/{env('GITHUB_RUN_ID')}"
    api(f"repos/{repo}/issues/{number}/comments", "POST",
        {"body": f"The nightly run is green again at {env('GITHUB_SHA')}: {run_url}"})
    api(f"repos/{repo}/issues/{number}", "PATCH",
        {"state": "closed", "state_reason": "completed"})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["run", "failure", "step", "success"])
    parser.add_argument("name", nargs="?")
    parser.add_argument("--command", action="append", default=[])
    parser.add_argument("--jobs", default="",
                        help="failure: triage only these comma-separated jobs")
    args = parser.parse_args(argv)
    if args.action == "run" and not args.name:
        parser.error("run needs a job name")
    return {"run": cmd_run, "failure": cmd_failure, "step": cmd_step,
            "success": cmd_success}[args.action](args) or 0


if __name__ == "__main__":
    sys.exit(main())
