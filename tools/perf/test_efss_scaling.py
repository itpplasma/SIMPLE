"""Focused parsing checks; run with python -m unittest discover -s tools/perf."""

import tempfile
from pathlib import Path
import unittest
from types import SimpleNamespace

import numpy as np

from efss_scaling import check_traced, phase_seconds, physics, run
from compare_efss_runs import compare


class MeasurementParsing(unittest.TestCase):
    def test_cpu_and_gpu_timers(self):
        self.assertEqual(
            phase_seconds(
                "INFO: Particle tracing completed      12.345 s (Total: 15.000 s)",
                False,
            ),
            12.345,
        )
        self.assertEqual(phase_seconds(" tracing        =     2.125000 s", True), 2.125)
        with self.assertRaises(ValueError):
            phase_seconds("STOP before tracing", False)

    def test_skipped_markers_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            np.savetxt(directory / "orbit_exit_code.dat", [[1, 0], [2, 1]])
            check_traced(directory, 2)
            with self.assertRaises(ValueError):
                check_traced(directory, 3)
            np.savetxt(directory / "orbit_exit_code.dat", [[1, 0], [2, 3]])
            with self.assertRaises(ValueError):
                check_traced(directory, 2)

    def test_fraction_semantics(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            # The passing/trapped columns are both normalized to resolved markers.
            np.savetxt(directory / "confined_fraction.dat", [[0, 0.25, 0.5, 8]])
            np.savetxt(directory / "unresolved_fraction.dat", [[0, 0.2, 10]])
            grid, confined, unresolved = physics(directory)
            np.testing.assert_array_equal(grid, [0])
            np.testing.assert_array_equal(confined, [0.75])
            np.testing.assert_array_equal(unresolved, [0.2])
            np.savetxt(directory / "confined_fraction.dat", [[0, np.nan, 0.5, 0]])
            with self.assertRaises(ValueError):
                physics(directory)

    def test_fraction_time_grids_must_match(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            np.savetxt(directory / "confined_fraction.dat", [[0, 0.5, 0.5, 10]])
            np.savetxt(directory / "unresolved_fraction.dat", [[1, 0, 10]])
            with self.assertRaises(ValueError):
                physics(directory)

    def test_existing_campaign_is_not_modified(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            metadata = directory / "metadata.json"
            metadata.write_text("original campaign")
            args = SimpleNamespace(
                executable=Path("unused"), wout=Path("unused"), output=directory
            )
            with self.assertRaises(FileExistsError):
                run(args)
            self.assertEqual(metadata.read_text(), "original campaign")


class BackendComparison(unittest.TestCase):
    def test_status_tolerance_and_start_identity(self):
        with tempfile.TemporaryDirectory() as name:
            cpu, gpu = Path(name) / "cpu", Path(name) / "gpu"
            starts = np.arange(100).reshape(20, 5)
            for directory in (cpu, gpu):
                directory.mkdir()
                np.savetxt(directory / "start.dat", starts)
                np.savetxt(
                    directory / "orbit_exit_code.dat",
                    np.column_stack((np.arange(1, 21), np.zeros(20))),
                )
                np.savetxt(directory / "confined_fraction.dat", [[0, 0.4, 0.4, 20]])
                np.savetxt(directory / "unresolved_fraction.dat", [[0, 0, 20]])
            self.assertTrue(compare(cpu, gpu)["passed"])
            np.savetxt(gpu / "confined_fraction.dat", [[0, 0.4, 0.2, 20]])
            self.assertFalse(compare(cpu, gpu)["passed"])
            np.savetxt(gpu / "start.dat", starts + 1)
            with self.assertRaises(AssertionError):
                compare(cpu, gpu)

    def test_unresolved_cannot_hide_disagreement(self):
        with tempfile.TemporaryDirectory() as name:
            cpu, gpu = Path(name) / "cpu", Path(name) / "gpu"
            for directory in (cpu, gpu):
                directory.mkdir()
                np.savetxt(directory / "start.dat", np.zeros((20, 5)))
                np.savetxt(
                    directory / "orbit_exit_code.dat",
                    np.column_stack((np.arange(1, 21), np.zeros(20))),
                )
                np.savetxt(directory / "confined_fraction.dat", [[0, 0.4, 0.4, 20]])
                np.savetxt(directory / "unresolved_fraction.dat", [[0, 0, 20]])
            np.savetxt(gpu / "unresolved_fraction.dat", [[0, 0.2, 20]])
            self.assertFalse(compare(cpu, gpu)["passed"])

    def test_one_marker_and_skipped_outputs(self):
        with tempfile.TemporaryDirectory() as name:
            cpu, gpu = Path(name) / "cpu", Path(name) / "gpu"
            for directory in (cpu, gpu):
                directory.mkdir()
                np.savetxt(directory / "start.dat", np.zeros((1, 5)))
                np.savetxt(directory / "orbit_exit_code.dat", [[1, 0]])
                np.savetxt(directory / "confined_fraction.dat", [[0, 0.5, 0.5, 1]])
                np.savetxt(directory / "unresolved_fraction.dat", [[0, 0, 1]])
            self.assertEqual(compare(cpu, gpu)["markers"], 1)
            np.savetxt(gpu / "orbit_exit_code.dat", [[1, 3]])
            with self.assertRaises(ValueError):
                compare(cpu, gpu)

    def test_backend_grid_lengths_must_match(self):
        with tempfile.TemporaryDirectory() as name:
            cpu, gpu = Path(name) / "cpu", Path(name) / "gpu"
            for directory in (cpu, gpu):
                directory.mkdir()
                np.savetxt(directory / "start.dat", np.zeros((1, 5)))
                np.savetxt(directory / "orbit_exit_code.dat", [[1, 0]])
            np.savetxt(cpu / "confined_fraction.dat", [[0, 0.5, 0.5, 1]])
            np.savetxt(cpu / "unresolved_fraction.dat", [[0, 0, 1]])
            np.savetxt(
                gpu / "confined_fraction.dat", [[0, 0.5, 0.5, 1], [1, 0.5, 0.5, 1]]
            )
            np.savetxt(gpu / "unresolved_fraction.dat", [[0, 0, 1], [1, 0, 1]])
            with self.assertRaises(ValueError):
                compare(cpu, gpu)


if __name__ == "__main__":
    unittest.main()
