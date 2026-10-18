#!/usr/bin/env python3
"""Export compact text sources from a completed public QA figure run."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from make_front import load_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--libneo-revision", required=True)
    args = parser.parse_args()
    results, confinement = load_results(args.workdir)
    exits = results["orbit_exit_code"]
    if np.any(exits >= 101):
        raise ValueError("Numerical failures in reference figure run")
    lost = exits == 1
    indices = np.flatnonzero(lost)
    rows = np.column_stack(
        [
            indices,
            results["zend"][lost, :3],
            results["xend_cart"][lost],
            results["times_lost"][lost],
            results["trap_par"][lost] >= 0,
        ]
    )
    args.output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output / "losses.csv",
        rows,
        delimiter=",",
        comments="",
        header="marker,s,theta_rad,phi_rad,x_cm,y_cm,z_cm,time_s,trapped",
        fmt=["%d"] + ["%.17g"] * 7 + ["%d"],
    )
    samples = np.unique(
        np.r_[
            np.arange(0, len(confinement), max(1, len(confinement) // 100)),
            len(confinement) - 1,
        ]
    )
    np.savetxt(
        args.output / "confinement.csv",
        confinement[samples, :3],
        delimiter=",",
        comments="",
        header="time_s,passing,trapped",
        fmt="%.17g",
    )
    metadata = {
        "format_version": 1,
        "markers": len(exits),
        "passing_markers": int(np.count_nonzero(results["trap_par"] < 0)),
        "losses": len(indices),
        "numerical_failures": 0,
        "plot_label": "public QA reference run",
        "reported_binary_version": results["version"],
        "reviewed_source_revision": args.source_revision,
        "libneo_revision": args.libneo_revision,
        "source_results_sha256": hashlib.sha256(
            (args.workdir / "results.nc").read_bytes()
        ).hexdigest(),
        "description": "Loss events and population counts for the public QA figure; "
        "confinement samples independently check the event-based curves.",
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
