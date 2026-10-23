"""Temporarily allow a leased main repair, retaining the exact protection state."""
import json
import os
from pathlib import Path
import sys
import urllib.request

URL = "https://gitlab.eufus.psnc.pl/api/v4/projects/itpplasma%2Fsimple/protected_branches/main"
SNAPSHOT = Path(os.environ["RUNNER_TEMP"]) / "eufus-main-protection.json"


def request(method="GET", payload=None):
    data = None if payload is None else json.dumps(payload).encode()
    req = urllib.request.Request(
        URL,
        data=data,
        method=method,
        headers={
            "PRIVATE-TOKEN": os.environ["EUFUS_GITLAB_TOKEN"],
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.load(response)


def main():
    if sys.argv[1] == "enable":
        original = request()
        # Save before mutation so the always() cleanup can recover a failed step.
        SNAPSHOT.write_text(json.dumps(original))
        SNAPSHOT.chmod(0o600)
        if not original["allow_force_push"]:
            changed = request("PATCH", {"allow_force_push": True})
            if changed != dict(original, allow_force_push=True):
                raise RuntimeError("EUFUS protection changed beyond the force flag")
        print("Saved main protection; force flag enabled for this dispatch")
    elif sys.argv[1] == "restore":
        if not SNAPSHOT.exists():
            print("No saved protection: no mutation was attempted")
            return
        original = json.loads(SNAPSHOT.read_text())
        current = request()
        if current["allow_force_push"] != original["allow_force_push"]:
            request("PATCH", {"allow_force_push": original["allow_force_push"]})
        if request() != original:
            raise RuntimeError("EUFUS main protection differs from its saved state")
        print("Verified exact restoration of EUFUS main protection")
    else:
        raise ValueError("Expected enable or restore")


if __name__ == "__main__":
    main()
