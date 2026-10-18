"""Archive existing completed marker CSVs; does not change numerical values."""
from pathlib import Path
import gzip
import json

ROOT = Path("concept/interface_candidates/results")

for suite, expected in (("main", 318), ("long", 29), ("sampling", 26)):
    manifest = ROOT / f"{suite}_manifest.json"
    if not manifest.exists():
        raise SystemExit(f"{suite} suite is absent; run study first")
    records = json.loads(manifest.read_text())
    if len(records) != expected:
        raise SystemExit(f"{suite} suite is incomplete ({len(records)}/{expected})")
    for rec in records:
        path = Path(rec["file"])
        if path.suffix == ".gz":
            continue
        archive = path.with_suffix(".csv.gz")
        archive.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        path.unlink()
        rec["file"] = str(archive)
    manifest.write_text(json.dumps(records, indent=2)+"\n")
print("Marker archives packed losslessly (gzip, deterministic timestamp).")
