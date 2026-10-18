#!/usr/bin/env python3
"""Check SIMPLE namelist keys, types, array bounds and basic value constraints.

Uses only the Python standard library and the committed input schema. This is
an input preflight; runtime field, file and cross-option checks remain in SIMPLE.
"""

from __future__ import annotations

import argparse
import difflib
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[de][+-]?\d+)?", re.IGNORECASE)
INTEGER = re.compile(r"[+-]?\d+")
TOKEN = re.compile(
    r"\s+|![^\n]*|'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"|[,=()*/:&$]|[^\s!,=()*/:&$'\"]+"
)


@dataclass
class Token:
    text: str
    line: int


class InputError(ValueError):
    pass


def tokenize(text: str) -> list[Token]:
    result = []
    position, line = 0, 1
    for match in TOKEN.finditer(text):
        if match.start() != position:
            raise InputError(f"line {line}: unterminated character value")
        word = match[0]
        if not word.isspace() and not word.startswith("!"):
            result.append(Token(word, line))
        line += word.count("\n")
        position = match.end()
    if position != len(text):
        raise InputError(f"line {line}: unterminated character value")
    return result


class Parser:
    def __init__(self, text: str):
        self.tokens = tokenize(text)
        self.i = 0

    def peek(self, offset: int = 0) -> str:
        pos = self.i + offset
        return self.tokens[pos].text.lower() if pos < len(self.tokens) else ""

    def take(self) -> Token:
        if not self.peek():
            raise InputError("unexpected end of input (missing '/'?)")
        token = self.tokens[self.i]
        self.i += 1
        return token

    def expect(self, word: str) -> None:
        token = self.take()
        if token.text.lower() != word:
            raise InputError(
                f"line {token.line}: expected {word!r}, got {token.text!r}"
            )

    def at_assignment(self) -> bool:
        return bool(re.fullmatch(r"[a-z_]\w*", self.peek())) and self.peek(1) in (
            "=",
            "(",
        )

    def ended(self) -> bool:
        return self.peek() == "/" or (
            self.peek() in ("&", "$") and self.peek(1) == "end"
        )

    def groups(self):
        while self.peek():
            marker = self.take()
            if marker.text not in ("&", "$"):
                raise InputError(
                    f"line {marker.line}: expected a namelist group, got {marker.text!r}"
                )
            group = self.take()
            assignments = []
            while not self.ended():
                if not self.peek():
                    raise InputError(
                        f"line {group.line}: unterminated &{group.text}; expected '/'"
                    )
                key = self.take()
                if not re.fullmatch(r"[a-z_]\w*", key.text, re.IGNORECASE):
                    raise InputError(
                        f"line {key.line}: expected a key, got {key.text!r}"
                    )
                indices = None
                if self.peek() == "(":
                    self.take()
                    indices = []
                    while self.peek() and self.peek() != ")":
                        indices.append(self.take().text)
                    self.expect(")")
                self.expect("=")
                supplied = []
                need_value = True
                while self.peek() and not self.ended() and not self.at_assignment():
                    if self.peek() == ",":
                        self.take()
                        if need_value:
                            supplied.append((None, 1))
                        need_value = True
                        continue
                    value = self.take()
                    count = 1
                    if self.peek() == "*":
                        self.take()
                        if not value.text.isdigit() or int(value.text) < 1:
                            raise InputError(
                                f"line {value.line}: repetition count must be positive"
                            )
                        count = int(value.text)
                        value = (
                            None
                            if (
                                self.peek() in (",", "")
                                or self.ended()
                                or self.at_assignment()
                            )
                            else self.take()
                        )
                    supplied.append((value.text if value else None, count))
                    need_value = False
                assignments.append((key, indices, supplied))
            if self.peek() == "/":
                self.take()
            else:
                self.take()
                self.expect("end")
            yield group, assignments


def scalar(value: str, type_name: str):
    if type_name == "integer" and INTEGER.fullmatch(value):
        return int(value)
    if type_name == "real" and NUMBER.fullmatch(value):
        number = float(re.sub(r"[dD]", "e", value))
        if math.isfinite(number):
            return number
    if type_name == "logical":
        # GNU namelist also accepts undotted T/F and TRUE/FALSE.
        if value.lower() in ("t", "true", ".true.", ".t."):
            return True
        if value.lower() in ("f", "false", ".false.", ".f."):
            return False
    if type_name == "character":
        if value[0] in "'\"" and value[-1] == value[0]:
            return value[1:-1].replace(value[0] * 2, value[0])
        # Undelimited character constants are valid namelist values.
        if re.fullmatch(r"[^\s,=()*/:&$'\"]+", value):
            return value
    raise InputError(f"expects {type_name}, got {value!r}")


def capacity(item: dict, indices: list[str] | None) -> int:
    extents = item["extent"]
    if indices is None:
        return math.prod(extents)
    if len(extents) != 1:
        raise InputError("subscripts require a rank-one array")
    text = "".join(indices)
    extent = extents[0]
    if INTEGER.fullmatch(text):
        index = int(text)
        if 1 <= index <= extent:
            # An element designator can start a sequence of array values.
            return extent - index + 1
        raise InputError(f"index {index} is outside 1:{extent}")
    parts = text.split(":")
    if 2 <= len(parts) <= 3 and all(not p or INTEGER.fullmatch(p) for p in parts):
        step = int(parts[2]) if len(parts) == 3 and parts[2] else 1
        if step == 0:
            raise InputError("array stride must be nonzero")
        # Fortran defaults to the declared bounds, regardless of stride sign.
        start = int(parts[0]) if parts[0] else 1
        end = int(parts[1]) if parts[1] else extent
        if not (1 <= start <= extent and 1 <= end <= extent):
            raise InputError(f"section {text} is outside 1:{extent}")
        return len(range(start, end + (1 if step > 0 else -1), step))
    raise InputError(f"unsupported array subscript {text!r}")


def check_rule(value, rule: dict) -> None:
    if "choices" in rule:
        choices = rule["choices"]
        candidate = value.lower() if isinstance(value, str) else value
        if candidate not in choices:
            raise InputError(f"must be one of {choices}, got {value!r}")
    for bound, operator in (("min", ">="), ("max", "<=")):
        if bound in rule and (
            value < rule[bound] if bound == "min" else value > rule[bound]
        ):
            raise InputError(f"must be {operator} {rule[bound]}")
    if rule.get("positive") and value <= 0:
        raise InputError("must be positive")
    if rule.get("nonnegative") and value < 0:
        raise InputError("must be nonnegative")


def check_file(path: Path, schema: dict) -> list[str]:
    errors = []
    seen_groups = set()
    try:
        for group_token, entries in Parser(path.read_text(encoding="utf-8")).groups():
            group = group_token.text.lower()
            seen_groups.add(group)
            if group not in schema["namelists"]:
                errors.append(f"{path}:{group_token.line}: unknown namelist '&{group}'")
                continue
            variables = {v["name"]: v for v in schema["namelists"][group]["variables"]}
            for key, indices, supplied in entries:
                name = key.text.lower()
                prefix = f"{path}:{key.line}: &{group} {name}"
                if name not in variables:
                    suggestion = difflib.get_close_matches(
                        name, variables, n=1, cutoff=0.55
                    )
                    hint = f"; did you mean '{suggestion[0]}'?" if suggestion else ""
                    errors.append(f"{prefix}: unknown key{hint}")
                    continue
                item = variables[name]
                try:
                    limit = capacity(item, indices)
                    count = sum(n for _, n in supplied)
                    if count > limit:
                        raise InputError(
                            f"accepts at most {limit} value(s), got {count}"
                        )
                    for raw, _ in supplied:
                        if raw is None:
                            continue  # Null values leave the compiled default unchanged.
                        value = scalar(raw, item["type"])
                        check_rule(value, item.get("validation", {}))
                except (InputError, OverflowError, ValueError) as error:
                    errors.append(f"{prefix}: {error}")
        if not ({"config", "poincare"} & seen_groups):
            errors.append(f"{path}: missing &config (or standalone &poincare) namelist")
    except (InputError, OSError, UnicodeError) as error:
        errors.append(f"{path}: {error}")
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs/input-reference.json",
    )
    args = parser.parse_args(argv)
    try:
        schema = json.loads(args.schema.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        parser.exit(1, f"Cannot read schema: {error}\n")
    errors = [error for path in args.inputs for error in check_file(path, schema)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    for path in args.inputs:
        print(f"OK {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
