"""Focused parsing checks; run with python -m unittest discover -s tools/perf."""
import tempfile
from pathlib import Path
import unittest

import numpy as np

from efss_scaling import phase_seconds, physics


class MeasurementParsing(unittest.TestCase):
    def test_cpu_and_gpu_timers(self):
        self.assertEqual(phase_seconds(
            'INFO: Particle tracing completed      12.345 s (Total: 15.000 s)', False), 12.345)
        self.assertEqual(phase_seconds(' tracing        =     2.125000 s', True), 2.125)
        with self.assertRaises(ValueError):
            phase_seconds('STOP before tracing', False)

    def test_fraction_semantics(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            # The passing/trapped columns are both normalized to resolved markers.
            np.savetxt(directory / 'confined_fraction.dat', [[0, .25, .5, 8]])
            np.savetxt(directory / 'unresolved_fraction.dat', [[0, .2, 10]])
            grid, confined, unresolved = physics(directory)
            np.testing.assert_array_equal(grid, [0])
            np.testing.assert_array_equal(confined, [.75])
            np.testing.assert_array_equal(unresolved, [.2])
            np.savetxt(directory / 'confined_fraction.dat', [[0, np.nan, .5, 0]])
            with self.assertRaises(ValueError):
                physics(directory)


if __name__ == '__main__':
    unittest.main()
