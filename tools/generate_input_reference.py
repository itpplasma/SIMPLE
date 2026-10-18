#!/usr/bin/env python3
"""Generate the SIMPLE namelist schema and input reference.

The source of truth is the Fortran ``namelist`` declaration.  This small
parser intentionally accepts the declaration subset used by SIMPLE and keeps
the generated files deterministic so that ``--check`` can be used in CI.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

RULES = {
    # Conservative checks from params::finalize_config and its callees.
    "gpu_mode": {"choices": ["off", "production", "benchmark"]},
    "gpu_backend": {"choices": ["openacc", "cuda-native", "auto"]},
    "gpu_method": {"choices": ["auto", "euler", "midpoint", "cash-karp", "dopri"]},
    "gpu_start_coordinates": {"choices": ["reference", "boozer"]},
    "wall_units": {"choices": ["m", "mm"]},
    "chart_boundary_kind": {"choices": ["auto", "lcfs", "wall", "domain"]},
    "macrostep_time_grid": {"choices": ["linear", "log"]},
    "axis_healing": {
        "choices": ["", "legacy", "legacy_adaptive", "powerlaw", "polyfit"]
    },
    "gpu_num_devices": {"min": 1},
    "gpu_t_block": {"positive": True},
    "gpu_loss_tau": {"positive": True},
    "gpu_maxloss": {"nonnegative": True},
    "gpu_min_timestep": {"nonnegative": True},
    "gpu_pilot_fraction": {"min": 0, "max": 0.25},
    "trace_time": {"positive": True},
    "relerr": {"positive": True},
    "canonical_ode_relerr": {"positive": True},
    "n_d": {"positive": True},
    "batch_size": {"min": 1},
    "face_al": {"positive": True},
    "ntestpart": {"min": 1},
    "nper": {"min": 1},
    "ntimstep": {"min": 2},
    "npoiper": {"min": 1},
    "npoiper2": {"min": 1},
    "notrace_passing": {"min": 0, "max": 1},
    "num_surf": {"min": 0, "max": 2},
    "startmode": {"min": 1, "max": 6},
    "integmode": {"choices": [-1, 0, 1, 2, 3, 4, 5, 6, 7, 15]},
    "orbit_model": {"choices": [0, 7]},
    "orbit_coord": {"choices": [0, 1]},
    "grid_density": {"min": 0, "max": 1},
    "cut_in_per": {"min": 0, "max": 1},
    "crossing_level": {"choices": [3]},
    "spectre_fo_entry_larmor": {"positive": True},
    "spectre_fo_tolerance": {"positive": True},
    "spectre_ncon_r": {"min": 1},
    "spectre_ncon_th": {"min": 1},
    "spectre_ncon_order": {"min": 1},
    "spectre_ncon_ode_max_steps": {"min": 1},
    "spectre_ncon_ode_relerr": {"positive": True},
    "canonical_grid_nr": {"min": 6, "max": 512},
    "canonical_grid_ntheta": {"min": 6, "max": 512},
    "canonical_grid_nphi": {"min": 6, "max": 512},
}


@dataclass
class Declaration:
    name: str
    type: str
    kind: str
    length: str
    shape: list[str]
    default: str
    source: str
    line: int
    comment: str


def split_comment(line: str) -> tuple[str, str]:
    quote = None
    i = 0
    while i < len(line):
        char = line[i]
        if char in "'\"":
            if quote == char:
                if i + 1 < len(line) and line[i + 1] == char:
                    i += 2
                    continue
                quote = None
            elif quote is None:
                quote = char
        elif char == "!" and quote is None:
            return line[:i], line[i + 1 :].strip().lstrip(">").strip()
        i += 1
    return line, ""


def split_top(text: str, delimiter: str = ",") -> list[str]:
    parts: list[str] = []
    start = 0
    depth = 0
    quote = None
    i = 0
    while i < len(text):
        char = text[i]
        if char in "'\"":
            if quote == char:
                if i + 1 < len(text) and text[i + 1] == char:
                    i += 2
                    continue
                quote = None
            elif quote is None:
                quote = char
        elif quote is None:
            if char in "([":
                depth += 1
            elif char in ")]":
                depth = max(0, depth - 1)
            elif char == delimiter and depth == 0:
                parts.append(text[start:i].strip())
                start = i + 1
        i += 1
    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


def statements(lines: list[str]) -> Iterable[tuple[int, str, str]]:
    """Join free-form continuations, retaining adjacent declaration comments."""
    code, comments, start = "", [], 0
    preceding: list[str] = []
    for number, raw in enumerate(lines, 1):
        part, comment = split_comment(raw)
        if not part.strip():
            if comment:
                (comments if code else preceding).append(comment)
            elif not code:
                preceding = []
            continue
        if not code:
            start = number
            comments = preceding
            preceding = []
        if comment:
            comments.append(comment)
        code += " " + part.strip().lstrip("&").rstrip().removesuffix("&")
        if part.rstrip().endswith("&"):
            continue
        yield start, code.strip(), " ".join(comments)
        code, comments = "", []
    if code:
        raise ValueError(f"line {start}: unfinished Fortran continuation")


def parse_decl_statement(
    code: str, source: str, line: int, comment: str
) -> list[Declaration]:
    match = re.match(
        r"(?is)^(integer|real|logical|character|complex|double\s+precision)\s*"
        r"(?:\(([^)]*)\)|\*(\d+))?\s*(?:,\s*[^:]*)?::\s*(.*)$",
        code,
    )
    if not match:
        return []
    base, selector, star, variables = match.groups()
    base = re.sub(r"\s+", " ", base.lower())
    kind, length = selector or star or "", ""
    if base == "character":
        length, kind = kind.removeprefix("len=").strip(), ""
    elif base == "double precision":
        base, kind = "real", "kind(1d0)"
    else:
        kind = kind.removeprefix("kind=").strip()
    dimension = re.search(r"(?is),\s*dimension\s*\(([^)]*)\)", code.split("::")[0])
    declaration_shape = split_top(dimension.group(1)) if dimension else []
    result = []
    for item in split_top(variables):
        var_match = re.fullmatch(
            r"(?is)([a-z_]\w*)\s*(?:\(([^)]*)\))?\s*(?:=\s*(.*))?", item
        )
        if not var_match:
            raise ValueError(f"{source}:{line}: unsupported declaration: {item}")
        name, shape, default = var_match.groups()
        result.append(
            Declaration(
                name.lower(),
                base,
                kind,
                length,
                split_top(shape) if shape else declaration_shape,
                (default or "").strip(),
                source,
                line,
                comment,
            )
        )
    return result


@dataclass
class Scope:
    declarations: dict[str, Declaration]
    uses: dict[str, tuple[str, str]]
    groups: list[tuple[str, int, list[str]]]
    source: str


def read_scopes(path: Path, source: str) -> dict[str, Scope]:
    """Index module/program specification parts, never unrelated local variables.

    Program-level literal assignments before the first CALL supply the defaults
    for spectre_poincare. No arbitrary later assignment is treated as a default.
    """
    scopes = {}
    scope = None
    specification = False
    in_type = False
    for line, code, comment in statements(
        path.read_text(encoding="utf-8").splitlines()
    ):
        lower = code.lower()
        begin = re.fullmatch(r"(?:module|program)\s+(\w+)", lower)
        if begin:
            scope = Scope({}, {}, [], source)
            scopes[begin[1]] = scope
            specification = True
            continue
        if re.match(r"end\s*(module|program)\b", lower):
            scope = None
        if scope is None:
            if lower.startswith("namelist") and not source.startswith("libneo/"):
                raise ValueError(
                    f"{source}:{line}: unsupported namelist scope; extend the extractor"
                )
            continue
        if lower == "contains" or re.match(r"call\s+", lower):
            specification = False
        if not specification:
            if lower.startswith("namelist") and not source.startswith("libneo/"):
                raise ValueError(
                    f"{source}:{line}: unsupported local namelist; extend the extractor"
                )
            continue
        if re.match(r"type\s*(?:,.*)?::\s*\w+", lower):
            in_type = True
        if re.match(r"end\s*type\b", lower):
            in_type = False
            continue
        if in_type:
            continue
        use = re.match(r"use\s+(\w+)\s*,\s*only\s*:\s*(.*)", lower)
        if use:
            for item in split_top(use[2]):
                names = [name.strip() for name in item.split("=>")]
                scope.uses[names[0]] = (use[1], names[-1])
        for decl in parse_decl_statement(code, source, line, comment):
            scope.declarations[decl.name] = decl
        group = re.fullmatch(r"namelist\s*/\s*(\w+)\s*/\s*(.*)", lower)
        if group:
            names = split_top(group[2])
            if not all(re.fullmatch(r"[a-z_]\w*", name) for name in names):
                raise ValueError(f"{source}:{line}: unsupported namelist declaration")
            scope.groups.append((group[1], line, names))
        elif lower.startswith("namelist"):
            raise ValueError(f"{source}:{line}: unsupported namelist declaration")
        assignment = re.fullmatch(r"([a-z_]\w*)\s*=\s*(.+)", code, re.IGNORECASE)
        if assignment and assignment[1].lower() in scope.declarations:
            decl = scope.declarations[assignment[1].lower()]
            if not decl.default:
                decl.default = assignment[2]
    return scopes


def build_schema(root: Path, libneo_source: Path | None = None) -> dict:
    libneo_source = libneo_source or root / "build/_deps/libneo-src"
    scopes = {}
    local_scopes = []
    for directory, prefix in (
        (root / "src", "src"),
        (root / "app", "app"),
        (libneo_source / "src", "libneo/src"),
    ):
        for path in sorted(directory.rglob("*")):
            if path.suffix.lower() not in (".f90", ".f"):
                continue
            # Fixed-form sources are not part of SIMPLE's current input
            # declaration graph.
            if path.suffix.lower() == ".f":
                continue
            parsed = read_scopes(path, f"{prefix}/{path.relative_to(directory)}")
            scopes.update(parsed)
            if prefix != "libneo/src":
                local_scopes.extend(parsed.values())

    def resolve(
        scope: Scope, name: str, visited: tuple = ()
    ) -> tuple[Scope, Declaration]:
        if (scope.source, name) in visited:
            raise ValueError(f"cyclic USE resolution for {name}")
        if name in scope.declarations:
            return scope, scope.declarations[name]
        if name in scope.uses:
            module, original = scope.uses[name]
            if module in scopes:
                return resolve(
                    scopes[module], original, visited + ((scope.source, name),)
                )
        raise ValueError(
            f"{scope.source}: cannot resolve declaration of {name}; "
            "configure SIMPLE first or pass --libneo-source PATH"
        )

    def constant(scope: Scope, text: str) -> str:
        if re.fullmatch(r"[a-z_]\w*", text, re.IGNORECASE):
            owner, decl = resolve(scope, text.lower())
            return constant(owner, decl.default)
        return text

    groups = {}
    for scope in local_scopes:
        for group, line, names in scope.groups:
            if group in groups:
                raise ValueError(
                    f"duplicate namelist &{group}; extend scope resolution"
                )
            entries = []
            for name in names:
                owner, declaration = resolve(scope, name)
                if not declaration.default:
                    raise ValueError(
                        f"{declaration.source}:{declaration.line}: "
                        f"no initial default found for {name}"
                    )
                shape = declaration.shape
                item = {
                    "name": name,
                    "type": declaration.type,
                    "kind": declaration.kind,
                    "length": declaration.length,
                    "shape": shape,
                    "extent": [int(constant(owner, size)) for size in shape],
                    "default": declaration.default,
                    "default_value": constant(owner, declaration.default),
                    "documentation": declaration.comment,
                    "source": declaration.source,
                    "source_line": declaration.line,
                }
                # Rules belong to &config, not to same-named app variables.
                if group == "config" and name in RULES:
                    item["validation"] = RULES[name]
                entries.append(item)
            groups[group] = {
                "source": scope.source,
                "source_line": line,
                "variables": entries,
            }
    return {
        "schema_version": 1,
        "generated_by": "tools/generate_input_reference.py",
        "namelists": groups,
    }


def render_rst(schema: dict) -> str:
    out = [
        "SIMPLE Input Reference",
        "======================",
        "",
        "Generated from the Fortran declarations; do not edit this page by hand.",
        "``simple.in`` requires ``&config`` and optionally contains ``&profiles``.",
        "``&poincare`` belongs to the separate ``spectre_poincare.in`` input.",
        "The readers are ``src/params.f90:read_config_namelist``,",
        "``src/profiles.f90:read_profiles`` and ``app/spectre_poincare.f90:read_input``.",
        "",
        "Validate an input from the repository root (Python standard library only)::",
        "",
        "   python tools/check_input.py simple.in examples/*.in",
        "",
        "The preflight reports unknown keys with suggestions, type errors, array",
        "bounds and selected basic value constraints from the JSON schema. It returns",
        "a nonzero status on errors. Runtime checks of files, fields and interacting",
        "options remain the responsibility of SIMPLE. This is not a full Fortran parser.",
        "It supports scalar and rank-one indexed/section assignments, repetition,",
        "null values, comments, quoted strings, and slash or END group terminators.",
        "",
        "Regenerate after configuring SIMPLE so that libneo sources are available::",
        "",
        "   make configure CONFIG=Fast",
        "   python tools/generate_input_reference.py",
        "   python tools/generate_input_reference.py --check",
        "",
        "For a custom libneo checkout, pass ``--libneo-source PATH`` to the generator.",
        "CTest checks both this page and the :download:`JSON schema <input-reference.json>` for staleness.",
        "Type selectors, shapes and initial defaults retain their source expressions;",
        "named constants are also resolved in JSON's ``default_value`` and ``extent``.",
        "A scalar initializer for an array applies to every element. Character lengths",
        "are shown with the type. Defaults precede runtime alias resolution and depend",
        "on the libneo sources selected for generation. An empty documentation cell",
        "means no adjacent declaration comment was found; see ``examples/simple_full.in``",
        "and ``DOC/config.md`` for additional usage, units and option interactions.",
        "",
    ]
    for group, data in schema["namelists"].items():
        out += [
            f"&{group}",
            "-" * (len(group) + 1),
            "",
            f"Source: ``{data['source']}:{data['source_line']}``.",
            "",
            ".. list-table::",
            "   :header-rows: 1",
            "   :widths: 15 10 10 10 20 35",
            "",
            "   * - Name",
            "     - Type",
            "     - Kind",
            "     - Shape",
            "     - Initial default",
            "     - Declaration documentation",
        ]
        for item in data["variables"]:
            shape = ", ".join(item["shape"]) if item["shape"] else "scalar"
            type_name = item["type"]
            if item["length"]:
                type_name += f"(len={item['length']})"
            doc = item["documentation"]
            for char in ("\\", "*", "`", "|"):
                doc = doc.replace(char, "\\" + char)
            out += [
                f"   * - ``{item['name']}``",
                f"     - ``{type_name}``",
                f"     - ``{item['kind'] or 'default'}``",
                f"     - ``{shape}``",
                f"     - ``{item['default']}``",
                f"     - {doc}",
            ]
        out.append("")
    return "\n".join(line.rstrip() for line in out).rstrip() + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--libneo-source", type=Path, help="libneo source used by the SIMPLE build"
    )
    parser.add_argument("--schema", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument(
        "--check",
        action="store_true",
        help="check committed outputs without modifying them",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    schema_path = (args.schema or root / "docs/input-reference.json").resolve()
    reference_path = (args.reference or root / "docs/input-reference.rst").resolve()
    try:
        schema = build_schema(root, args.libneo_source)
    except (ValueError, OSError) as error:
        parser.exit(1, f"{error}\n")
    schema_text = json.dumps(schema, indent=2, sort_keys=False) + "\n"
    reference_text = render_rst(schema)
    if args.check:
        stale = []
        if (
            not schema_path.exists()
            or schema_path.read_text(encoding="utf-8") != schema_text
        ):
            stale.append(str(schema_path))
        if (
            not reference_path.exists()
            or reference_path.read_text(encoding="utf-8") != reference_text
        ):
            stale.append(str(reference_path))
        if stale:
            print("Input reference is stale: " + ", ".join(stale), file=sys.stderr)
            print(
                "Regenerate with: python tools/generate_input_reference.py",
                file=sys.stderr,
            )
            return 1
        return 0
    schema_path.parent.mkdir(parents=True, exist_ok=True)
    schema_path.write_text(schema_text, encoding="utf-8")
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    reference_path.write_text(reference_text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
