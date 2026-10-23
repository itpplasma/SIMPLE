"""Exercise golden executable selection without compiling or downloading data."""

from pathlib import Path
import os
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
HARNESS = ROOT / "test/golden_record/golden_record.sh"
LIBNEO_REF = "924ba7735366d95c2c6855bf00bd0e46b8c05cf2"


def executable(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\nset -eu\n" + body)
    path.chmod(0o755)


def main():
    with tempfile.TemporaryDirectory(prefix="golden-build-") as tmp:
        base = Path(tmp)
        reference, current = base / "reference", base / "current"
        alternate = base / "alternate build"
        data = base / "golden/test_data"
        data.mkdir(parents=True)
        for name in ("wout.nc", "wout_ncsx.nc", "coils.c09r00.simple"):
            (data / name).touch()
        # Only rev-parse reads this repository; the fixture never changes it.
        git_dir = subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "--absolute-git-dir"], text=True
        ).strip()
        reference.mkdir()
        (reference / ".git").write_text(f"gitdir: {git_dir}\n")
        for project in (reference, current):
            case = project / "test/golden_record/albert_coils"
            case.mkdir(parents=True)
            (case / "simple.in").write_text("&config\ntrace_time = 1d-4\n/\n")
        for build in (reference / "build", alternate):
            executable(
                build / "simple.x",
                """
cat > times_lost.dat <<'EOF'
1 1e-4 0.5 0.3 0.1 0.3 0 0 1 0
2 1e-4 0.5 0.3 0.1 0.3 0 0 1 0
EOF
printf '1 0 1e-4 0 0\n2 0 1e-4 0 0\n' > orbit_exit_code.dat
printf '0 1\n1e-4 1\n' > confined_fraction.dat
""",
            )
            for program, output in (
                (
                    "field_can/test_field_can_albert_coils_diagnostic.x",
                    "albert_coils_diagnostic.dat",
                ),
                ("field_can/test_albert_intermediate.x", "albert_intermediate.dat"),
                ("test_sympl_tok.x", "selected_tokamak.dat"),
            ):
                body = f"printf '1 2\\n' > {output}\n"
                if program == "test_sympl_tok.x":
                    body += (
                        "for method in euler1 euler2 midpoint gauss4; do\n"
                        "  printf '1 2\\n' > tokamak_testfield_$method.dat\n"
                        "done\n"
                    )
                executable(build / "test/tests" / program, body)
        (reference / "build/.golden-libneo-ref").write_text(LIBNEO_REF + "\n")
        executable(current / "build/simple.x", "exit 97\n")
        env = os.environ.copy()
        env.update(
            GOLDEN_RECORD_BASE_DIR=str(base / "golden"),
            GOLDEN_RECORD_REF_DIR=str(reference),
            GOLDEN_RECORD_CUR_DIR=str(current),
            GOLDEN_RECORD_CUR_BUILD_DIR=str(alternate),
            GOLDEN_LIBNEO_REF=LIBNEO_REF,
            GOLDEN_RECORD_OMP_THREADS="1",
        )
        result = subprocess.run(
            [str(HARNESS), "main", "albert_coils"],
            env=env,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        assert result.returncode == 0, result.stdout
        candidate_runs = [
            p for p in (base / "golden/runs").iterdir() if p.name != "run_main"
        ]
        assert len(candidate_runs) == 1, result.stdout
        run = candidate_runs[0]
        for output in ("albert_coils_diagnostic.dat", "albert_intermediate.dat"):
            assert (run / "albert_coils" / output).read_text() == "1 2\n"
        assert (run / "sympl_tokamak/selected_tokamak.dat").read_text() == "1 2\n"
        print("Alternate golden build selected for SIMPLE and diagnostics.")


if __name__ == "__main__":
    main()
