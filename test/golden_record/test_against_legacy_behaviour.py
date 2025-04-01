import sys
import numpy as np

OLD_FILE = sys.argv[1]
NEW_FILE = sys.argv[2]

old_lines = []
new_lines = []

old_data = np.loadtxt(OLD_FILE)
new_data = np.loadtxt(NEW_FILE)

non_match_count = 0
for old, new in zip(old_data, new_data):
    for old_d, new_d in zip(old, new):
        if not np.isclose(old_d, new_d):
            non_match_count += 1
print(
    f"Found {non_match_count} non-matching entries comparing {OLD_FILE} and {NEW_FILE}."
)
