"""Wrap emitted statements without changing their tokens or arithmetic."""

from pathlib import Path
import sys
import textwrap

for path in Path(sys.argv[1]).glob("fo_poly_*.inc"):
    lines = []
    logical = []
    pending = ""
    for line in path.read_text().splitlines():
        pending += (" " if pending else "") + line.strip().removesuffix("&").lstrip(
            "&"
        ).strip()
        if not line.strip().endswith("&"):
            logical.append(pending)
            pending = ""
    if pending:
        raise ValueError("unterminated emitted continuation")
    for line in logical:
        if line.startswith("!"):
            lines.append(line)
            continue
        wrapped = textwrap.wrap(
            line, width=80, break_long_words=False, break_on_hyphens=False
        ) or [""]
        lines.extend(
            text + (" &" if i < len(wrapped) - 1 else "")
            for i, text in enumerate(wrapped)
        )
    path.write_text("\n".join(lines) + "\n")
