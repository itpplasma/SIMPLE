#!/usr/bin/env python3
"""Check native Albert diagnostics without a Python plotting dependency."""

from __future__ import annotations

from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import zlib


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_png(path: Path) -> tuple[int, int, list[bytearray]]:
    """Decode the native backend's noninterlaced, eight-bit RGB PNG output."""
    source = path.read_bytes()
    require(source[:8] == b"\x89PNG\r\n\x1a\n", f"{path.name}: invalid PNG signature")
    position = 8
    compressed = bytearray()
    header = None
    ended = False
    while position < len(source):
        length = struct.unpack_from(">I", source, position)[0]
        kind = source[position + 4:position + 8]
        chunk = source[position + 8:position + 8 + length]
        crc = struct.unpack_from(">I", source, position + 8 + length)[0]
        require(zlib.crc32(kind + chunk) == crc, f"{path.name}: corrupt PNG chunk")
        position += length + 12
        if kind == b"IHDR":
            header = struct.unpack(">IIBBBBB", chunk)
        elif kind == b"IDAT":
            compressed.extend(chunk)
        elif kind == b"IEND":
            ended = True
            break
    require(ended and header is not None, f"{path.name}: incomplete PNG")
    width, height, depth, color, compression, filtering, interlace = header
    require((depth, color, compression, filtering, interlace) == (8, 2, 0, 0, 0),
            f"{path.name}: expected native eight-bit RGB PNG")
    raw = zlib.decompress(compressed)
    stride = 3 * width
    require(len(raw) == height * (stride + 1), f"{path.name}: invalid image data size")
    rows = []
    previous = bytearray(stride)
    for y in range(height):
        offset = y * (stride + 1)
        filter_type = raw[offset]
        require(filter_type <= 4, f"{path.name}: invalid PNG filter")
        row = bytearray(raw[offset + 1:offset + stride + 1])
        if filter_type:
            for x in range(stride):
                left = row[x - 3] if x >= 3 else 0
                above = previous[x]
                upper_left = previous[x - 3] if x >= 3 else 0
                if filter_type == 1:
                    prediction = left
                elif filter_type == 2:
                    prediction = above
                elif filter_type == 3:
                    prediction = (left + above) // 2
                else:
                    p = left + above - upper_left
                    candidates = (left, above, upper_left)
                    prediction = min(candidates, key=lambda value: abs(p - value))
                row[x] = (row[x] + prediction) & 255
        rows.append(row)
        previous = row
    return width, height, rows


def main() -> None:
    require(len(sys.argv) == 2, "usage: test_diag_albert.py <Fortran test executable>")
    executable = Path(sys.argv[1]).resolve(strict=True)
    # CTest supplies the executable; this wrapper runs only as a project test.
    with tempfile.TemporaryDirectory(prefix="diag-albert-", dir=Path.cwd()) as temporary:
        directory = Path(temporary)
        subprocess.run([str(executable)], cwd=directory, check=True, timeout=120)
        expected_names = {
            f"albert_{field}_{slice_name}_contour.png"
            for field in ("Aph", "hth", "hph", "Bmod")
            for slice_name in ("inner", "middle", "outer")
        }
        actual_names = {path.name for path in directory.glob("albert_*_contour.png")}
        require(actual_names == expected_names, "Albert diagnostic output filenames changed")
        require(not list(directory.glob("*.py")), "native diagnostics wrote Python sidecars")
        for name in sorted(expected_names):
            actual = read_png(directory / name)
            expected = read_png(directory / ("expected_" + name))
            require(actual == expected,
                    f"{name}: field, radial slice, coordinates, or title changed")
        for name in ("albert_rectangular.png", "albert_square.png"):
            read_png(directory / name)
    print("Albert native PNG checks passed for all twelve field/slice diagnostics.")


if __name__ == "__main__":
    main()
