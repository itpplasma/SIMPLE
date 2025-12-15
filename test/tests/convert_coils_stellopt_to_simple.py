#!/usr/bin/env python3
"""Convert STELLOPT coils format to SIMPLE format."""

import sys

if len(sys.argv) < 3:
    print("Usage: convert_coils_stellopt_to_simple.py <input.coils> <output.simple>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as f:
    lines = f.readlines()

coords = []
in_filament = False

for line in lines:
    line = line.strip()
    if 'begin filament' in line.lower():
        in_filament = True
        continue
    if not in_filament or not line or line.startswith('#'):
        continue
    parts = line.split()
    if len(parts) >= 4:
        try:
            x = float(parts[0])
            y = float(parts[1])
            z = float(parts[2])
            current = float(parts[3])
            coords.append((x, y, z, current))
        except ValueError:
            pass

with open(output_file, 'w') as f:
    f.write(str(len(coords)) + '\n')
    for x, y, z, current in coords:
        f.write(f'{x} {y} {z} {current}\n')

print(f"Converted {len(coords)} coil segments from {input_file} to {output_file}")
