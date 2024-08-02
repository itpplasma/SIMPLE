
import numpy as np

OLD_FILE = 'times_lost_old.dat'
NEW_FILE = 'times_lost_new.dat'

old_lines = []
new_lines = []

old_data = np.loadtxt(OLD_FILE)
new_data = np.loadtxt(NEW_FILE)

i = 0
for old, new in zip(old_data, new_data):
    for old_d, new_d in zip(old, new):
        if not np.isclose(old_d, new_d):
            i += 1
print(f"Found {i} non-matching entries when comparing times_lost.dat.")
