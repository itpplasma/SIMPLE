"""Tests of the pysimple API against simple.x and physical invariants."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pytest

pysimple = pytest.importorskip("pysimple")

SIMPLE_X = Path(__file__).resolve().parents[2] / "build" / "simple.x"


def _write_config(path: Path, vmec_file: str, **extra: str) -> Path:
    settings = {
        "netcdffile": f"'{Path(vmec_file).resolve()}'",
        "multharm": "3",
        "ntestpart": "8",
        "trace_time": "1d-4",
        "sbeg": "0.4",
        "contr_pp": "-1d10",
        "npoiper2": "64",
        "isw_field_type": "2",
        "integmode": "3",
        "deterministic": ".True.",
        **extra,
    }
    body = "\n".join(f"{key} = {value}" for key, value in settings.items())
    config = path / "simple.in"
    config.write_text(f"&config\n{body}\n/\n")
    return config


class TestConfiguration:
    def test_config_file_reproduces_simple_x(self, vmec_file, tmp_path, monkeypatch):
        if not SIMPLE_X.exists():
            pytest.skip("simple.x is not built")
        config = _write_config(tmp_path, vmec_file)
        subprocess.run([str(SIMPLE_X), str(config)], cwd=tmp_path, check=True,
                       capture_output=True)
        reference_start = np.loadtxt(tmp_path / "start.dat").T
        reference_loss = np.loadtxt(tmp_path / "times_lost.dat")[:, 1]

        python_dir = tmp_path / "python"
        python_dir.mkdir()
        monkeypatch.chdir(python_dir)
        pysimple.init(config=config)
        start = pysimple.sample()
        result = pysimple.trace_parallel(start)

        np.testing.assert_allclose(start, reference_start, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(result["loss_times"], reference_loss,
                                   rtol=1e-10)

    def test_settings_override_config_file(self, vmec_file, tmp_path):
        config = _write_config(tmp_path, vmec_file)
        pysimple.init(config=config, trace_time=2e-4, sbeg=[0.3])

        assert pysimple.params.trace_time == 2e-4
        assert pysimple.params.ntestpart == 8
        assert pysimple.params.sbeg[0] == 0.3

    def test_init_starts_from_defaults(self, vmec_file):
        pysimple.init(vmec_file, multharm=3, facE_al=700.0)
        assert pysimple.params.face_al == 700.0

        pysimple.init(vmec_file, multharm=3)
        assert pysimple.params.face_al == 1.0

    def test_unknown_setting_is_rejected(self, vmec_file):
        with pytest.raises(ValueError, match="no_such_setting"):
            pysimple.init(vmec_file, no_such_setting=1)

    def test_missing_config_file_is_rejected(self, tmp_path):
        with pytest.raises(ValueError, match="missing.in"):
            pysimple.init(config=tmp_path / "missing.in")


class TestChartmap:
    def test_chartmap_selects_boozer_field(self):
        chartmap = os.environ.get("SIMPLE_TEST_CHARTMAP")
        if not chartmap:
            pytest.skip("SIMPLE_TEST_CHARTMAP not set")
        pysimple.init(chartmap, deterministic=True, npoiper2=64, trace_time=1e-4)

        result = pysimple.compute_canonical_frequencies(
            np.array([0.4, 0.7, 0.1, 1.0, 0.1]), n_periods=1, max_steps=10
        )
        # A boundary stop is a physical loss, never an integrator error.
        assert result["status"] in {pysimple.FREQ_SUCCESS, pysimple.FREQ_ORBIT_LOST,
                                    pysimple.FREQ_MAX_STEPS}

    def test_chartmap_invariants_use_lcfs_flux(self, tmp_path):
        source = os.environ.get("SIMPLE_TEST_CHARTMAP")
        if not source:
            pytest.skip("SIMPLE_TEST_CHARTMAP not set")
        chartmap = tmp_path / "extended_chartmap.nc"
        shutil.copy2(source, chartmap)
        rho_lcfs = 0.8
        with nc.Dataset(chartmap, "r+") as dataset:
            dataset.setncattr("rho_lcfs", rho_lcfs)
            s = np.asarray(dataset["s"][:])
            a_phi = np.asarray(dataset["A_phi"][:])
        expected_edge_flux = (a_phi[-1] - a_phi[0]) / (s[-1] - s[0]) * rho_lcfs**2

        pysimple.init(chartmap, deterministic=True, npoiper2=64, trace_time=1e-4)
        invariant = pysimple.invariants_from_state(np.array([0.4, 0.7, 0.1, 1.0, 0.1]))

        assert invariant["edge_flux_native"] == pytest.approx(expected_edge_flux,
                                                              rel=1e-10)


class TestTracing:
    def test_volume_samples_lie_in_shell(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, multharm=3)
        start = pysimple.sample_volume(10, s_inner=0.2, s_outer=0.8)

        assert np.all((start[0] >= 0.2) & (start[0] <= 0.8))
        assert np.all(np.abs(start[4]) <= 1.0)

    def test_unknown_integrator_is_rejected(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, multharm=3)
        start = pysimple.sample_surface(1, s=0.4)
        with pytest.raises(ValueError, match="does_not_exist"):
            pysimple.trace_parallel(start, integrator="does_not_exist")

    def test_skipped_passing_particle_keeps_sentinel(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-4, multharm=3)
        particle = np.array([[0.6], [0.0], [0.25 * np.pi], [1.0], [0.8]])

        traced = pysimple.trace_parallel(particle, integrator="midpoint")
        classified = pysimple.classify_parallel(particle, integrator="midpoint")

        for result in (traced, classified):
            assert result["loss_times"][0] == -1.0
            np.testing.assert_allclose(result["final_positions"][:, 0], particle[:, 0])

    def test_single_orbit_matches_batch_loss_time(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, trace_time=1e-2, multharm=3)
        particle = np.array([0.6, 0.0, 0.25 * np.pi, 1.0, 0.2])

        single = pysimple.trace_orbit(particle, integrator="midpoint",
                                      return_trajectory=True)
        batch = pysimple.trace_parallel(particle.reshape(5, 1), integrator="midpoint")

        assert single["loss_time"] == pytest.approx(batch["loss_times"][0])
        assert 0.0 < single["loss_time"] < pysimple.params.trace_time
        last = np.flatnonzero(single["times"] > 0.0)[-1]
        assert np.isnan(single["times"][last + 1:]).all()


class TestOrbitDiagnostics:
    def test_frequency_arguments_are_validated(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, multharm=3)
        particle = np.array([0.4, 0.7, 0.1, 1.0, 0.9])

        with pytest.raises(ValueError, match="n_periods"):
            pysimple.compute_canonical_frequencies(particle, n_periods=0)
        with pytest.raises(ValueError, match="shape"):
            pysimple.compute_canonical_frequencies(particle[:4])

    def test_trace_to_tip_stops_at_turning_point(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, npoiper2=1024, multharm=3)
        result = pysimple.trace_to_cut(np.array([0.4, 0.7, 0.1, 1.0, 0.1]),
                                       cut="tip", max_events=10, integrator="midpoint")

        assert result["status"] == 0
        assert result["cut_type"] == "tip"
        assert abs(result["state"][4]) < 1e-6

    def test_trace_to_cut_arguments_are_validated(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, multharm=3)
        particle = np.array([0.4, 0.7, 0.1, 1.0, 0.9])
        with pytest.raises(ValueError, match="Unknown cut"):
            pysimple.trace_to_cut(particle, cut="poloidal")
        with pytest.raises(ValueError, match="symplectic"):
            pysimple.trace_to_cut(particle, integrator="rk45")

    def test_states_from_invariants_reproduce_invariants(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, isw_field_type=-1)
        particle = np.array([0.2, 0.7, 0.2, 0.8, 0.6])
        invariant = pysimple.invariants_from_state(particle)

        assert invariant["status"] == pysimple.INVARIANT_SUCCESS
        assert invariant["h0"] == pytest.approx(particle[3] ** 2)

        starts = pysimple.states_from_invariants(
            invariant["h0"], invariant["j_perp"], invariant["p_phi"]
        )
        assert starts["status"] == pysimple.INVARIANT_SUCCESS
        assert starts["n_solutions"] >= 1
        assert set(starts["sigma"]).issubset({-1, 1})
        for state in starts["states"]:
            reconstructed = pysimple.invariants_from_state(state)
            assert reconstructed["h0"] == pytest.approx(invariant["h0"])
            assert reconstructed["j_perp"] == pytest.approx(invariant["j_perp"],
                                                            rel=1e-10)
            assert reconstructed["p_phi"] == pytest.approx(invariant["p_phi"],
                                                           rel=1e-9, abs=1e-9)

    def test_potato_invariants_are_normalized(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, isw_field_type=-1)
        result = pysimple.states_from_potato_invariants(
            1.0, 2.0, 12.0, psi_axis=10.0, psi_edge=20.0, v0_ratio=0.5
        )
        np.testing.assert_allclose(result["simple_invariants"], [0.25, 0.5, 2.0])

        reversed_flux = pysimple.states_from_potato_invariants(
            1.0, 2.0, 8.0, psi_axis=10.0, psi_edge=0.0
        )
        np.testing.assert_allclose(reversed_flux["simple_invariants"], [1.0, 2.0, 2.0])

    def test_invariant_arguments_are_validated(self, vmec_file):
        pysimple.init(vmec_file, deterministic=True, isw_field_type=-1)
        with pytest.raises(ValueError, match="shape"):
            pysimple.invariants_from_state(np.zeros(4))
        with pytest.raises(ValueError, match="h0"):
            pysimple.states_from_invariants(0.0, 0.0, 0.0)
        with pytest.raises(ValueError, match="v0_ratio"):
            pysimple.states_from_potato_invariants(
                1.0, 1.0, 1.0, psi_axis=0.0, psi_edge=1.0, v0_ratio=0.0
            )
