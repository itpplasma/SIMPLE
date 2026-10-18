#!/usr/bin/env python3
"""Focused regression tests for source extraction and namelist preflight."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from check_input import Parser, check_file, scalar
from generate_input_reference import (
    build_schema,
    parse_decl_statement,
    split_comment,
    split_top,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads((ROOT / "docs/input-reference.json").read_text())


class InputValidationTests(unittest.TestCase):
    def check_text(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "simple.in"
            path.write_text(text)
            return check_file(path, SCHEMA)

    def test_repository_inputs(self):
        paths = [ROOT / "simple.in", *sorted((ROOT / "examples").rglob("*.in"))]
        for path in paths:
            with self.subTest(path=path.relative_to(ROOT)):
                self.assertEqual(check_file(path, SCHEMA), [])
        print(f"Validated {len(paths)} repository inputs")

    def test_valid_syntax(self):
        for text in (
            "&CONFIG NPER=10 /",
            "&config startmode=7, fusion_grid=128,32,32 /",
            "&config field_input='it''s!/a=b.nc', debug=F /",
            "&config sbeg(2:4)=3*0.2, sbeg(1024)=0.1 /",
            "&config sbeg(4:2:-1)=0.1 0.2 0.3 /",
            "&config sbeg(:1:-1)=0.1, sbeg(1024::-1)=0.2 /",
            "&config sbeg=,0.2,2*,0.3, nper=, /",
            "&config nper=10 nper=20 /",  # Fortran permits later assignments.
            "&config\n debug=.TRUE. ! slash / in a comment\n /",
            "$config nper=10 $end",
            "&config &end\n&profiles Te_coef(1:10)=1d0,9*0d0 /",
            "&poincare nturns=200 spectre_file='spectre.h5' /",
        ):
            with self.subTest(text=text):
                self.assertEqual(self.check_text(text), [])

    def test_invalid_syntax_and_values(self):
        for text, diagnostic in (
            ("&config nperr=2 /", "did you mean 'nper'"),
            ("&config nper=1.5 /", "expects integer"),
            ("&config debug=banana /", "expects logical"),
            ("&config relerr=1e999 /", "expects real"),
            ("&config trace_time=-1 /", "must be positive"),
            ("&config facE_al=0 /", "must be positive"),
            ("&config ntimstep=1 /", "must be >= 2"),
            ("&config startmode=8 /", "must be <= 7"),
            ("&config fusion_grid=128,32,32,16 /", "at most 3"),
            ("&config gpu_mode='typo' /", "must be one of"),
            ("&config nper=1,2 /", "at most 1"),
            ("&config sbeg(1025)=0.2 /", "outside 1:1024"),
            ("&config sbeg(1:2)=3*0.2 /", "at most 2"),
            ("&config sbeg(::-1)=0.1 /", "at most 0"),
            ("&config sbeg(:1:-1)=0.1,0.2 /", "at most 1"),
            ("&config sbeg(1023::-1)=0.1 /", "at most 0"),
            ("&config sbeg=999999999999*0.2 /", "at most 1024"),
            ("&config sbeg(1:2:0)=0.2 /", "stride"),
            ("&config sbeg=0*0.2 /", "repetition count"),
            ("&config nper=2", "unterminated &config"),
            ("&config field_input='unterminated /", "unterminated character"),
            ("&config nper 2 /", "expected '='"),
            ("&config / trailing", "expected a namelist"),
            ("&typo /", "unknown namelist"),
            ("", "missing &config"),
            ("&profiles /", "missing &config"),
        ):
            with self.subTest(text=text):
                self.assertTrue(
                    any(diagnostic in e for e in self.check_text(text)),
                    self.check_text(text),
                )

    def test_assignment_line_numbers(self):
        errors = self.check_text("&config\n nper=10\n debug=banana\n /")
        self.assertIn(":3: &config debug:", errors[0])

    def test_cli_status(self):
        command = [sys.executable, str(ROOT / "tools/check_input.py")]
        good = subprocess.run(
            command + [str(ROOT / "simple.in")],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(good.returncode, 0, good.stderr)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.in"
            path.write_text("&config nperr=2 nper='wrong' trace_time=-1 /\n")
            bad = subprocess.run(
                command + [str(path)], check=False, capture_output=True, text=True
            )
            self.assertEqual(bad.returncode, 1)
            self.assertIn("did you mean 'nper'", bad.stderr)
            self.assertIn("expects integer", bad.stderr)
            self.assertIn("must be positive", bad.stderr)


class ExtractionTests(unittest.TestCase):
    def test_declaration_metadata(self):
        decl = parse_decl_statement(
            "real(dp), dimension(3) :: a=0d0, b=1d0", "fixture.f90", 7, "units"
        )[1]
        self.assertEqual(
            (decl.name, decl.type, decl.kind, decl.shape, decl.default, decl.comment),
            ("b", "real", "dp", ["3"], "1d0", "units"),
        )
        char = parse_decl_statement("character*1000 :: file=''", "fixture.f90", 1, "")[
            0
        ]
        self.assertEqual((char.type, char.length), ("character", "1000"))
        self.assertEqual(split_comment("character(*) :: a='it''s!ok' ! doc")[1], "doc")

    def fixture(self, directory):
        root = Path(directory)
        (root / "src").mkdir()
        (root / "docs").mkdir()
        (root / "src/input.f90").write_text("""module input
 use imported, only: alias => actual
 implicit none
 integer, parameter :: n=3
 ! Grid weights
 real(8) :: weights(n)=0d0
 namelist /config/ weights, alias
contains
 subroutine read_input
 read(10,nml=config)
 end subroutine
end module
""")
        (root / "src/imported.f90").write_text("""module imported
 logical :: actual=.true. ! enable output
end module
""")
        return root

    def test_imports_comments_and_shapes(self):
        with tempfile.TemporaryDirectory() as directory:
            schema = build_schema(self.fixture(directory))
        weights, alias = schema["namelists"]["config"]["variables"]
        self.assertEqual(weights["extent"], [3])
        self.assertEqual(weights["documentation"], "Grid weights")
        self.assertEqual(
            (alias["name"], alias["type"], alias["default"], alias["documentation"]),
            ("alias", "logical", ".true.", "enable output"),
        )

    def test_unresolved_declaration_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = self.fixture(directory)
            (root / "src/imported.f90").unlink()
            with self.assertRaisesRegex(
                ValueError, "cannot resolve declaration of alias"
            ):
                build_schema(root)
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tools/generate_input_reference.py"),
                    "--root",
                    str(root),
                    "--check",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("cannot resolve declaration of alias", result.stderr)
            self.assertIn(
                "Regenerate with: python tools/generate_input_reference.py",
                result.stderr,
            )

    def test_source_change_makes_reference_stale(self):
        command = [sys.executable, str(ROOT / "tools/generate_input_reference.py")]
        with tempfile.TemporaryDirectory() as directory:
            root = self.fixture(directory)
            result = subprocess.run(
                command + ["--root", str(root)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            path = root / "src/input.f90"
            path.write_text(
                path.read_text()
                .replace("weights, alias", "weights, alias, new_key")
                .replace(
                    " namelist /config/", " integer :: new_key=4\n namelist /config/"
                )
            )
            result = subprocess.run(
                command + ["--root", str(root), "--check"],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("Input reference is stale", result.stderr)
            self.assertIn("python tools/generate_input_reference.py", result.stderr)

    def test_defaults_match_compiled_reference_input(self):
        # test_config_reference checks this independently maintained input
        # against the compiled Fortran defaults, including imported libneo keys.
        for group, entries in Parser(
            (ROOT / "examples/simple_full.in").read_text()
        ).groups():
            documented = {key.text.lower(): values for key, _, values in entries}
            for item in SCHEMA["namelists"][group.text.lower()]["variables"]:
                with self.subTest(group=group.text, key=item["name"]):
                    text = item["default_value"]
                    if text.startswith("["):
                        defaults = split_top(text[1:-1])
                    elif text.startswith("(/"):
                        defaults = split_top(text[2:-2])
                    else:
                        defaults = [text] * (item["extent"][0] if item["extent"] else 1)
                    # Kind suffixes occur in source initializers, not namelist data.
                    defaults = [
                        scalar(re.sub(r"_(dp|8)$", "", v), item["type"])
                        for v in defaults
                    ]
                    supplied = [
                        scalar(v, item["type"])
                        for v, n in documented[item["name"]]
                        for _ in range(n)
                    ]
                    self.assertEqual(defaults, supplied)

    def test_scope_specific_defaults(self):
        poincare = {v["name"]: v for v in SCHEMA["namelists"]["poincare"]["variables"]}
        config = {v["name"]: v for v in SCHEMA["namelists"]["config"]["variables"]}
        self.assertEqual(poincare["nturns"]["default"], "200")
        self.assertEqual(config["nturns"]["default"], "8")
        self.assertEqual(config["integmode"]["default_value"], "1")
        self.assertTrue(config["netcdffile"]["source"].startswith("libneo/"))


if __name__ == "__main__":
    unittest.main()
