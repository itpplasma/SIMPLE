"""Behavioural tests of the shipped read-only plotting commands."""

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from pysimple import cli


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_wilson_reference_and_boundaries():
    # Independently tabulated 95% Wilson limits for 5 losses out of 20.
    low, high = cli.wilson_interval(5, 20)
    assert low == pytest.approx(0.11186170140766569)
    assert high == pytest.approx(0.468700877618744)
    low, high = cli.wilson_interval([0, 10, 0], [10, 10, 0])
    assert low[0] == pytest.approx(0, abs=1e-15)
    assert high[0] == pytest.approx(0.2775327998628892)
    assert low[1] == pytest.approx(1 - high[0])
    assert high[1] == pytest.approx(1)
    assert np.isnan(low[2]) and np.isnan(high[2])
    narrow = cli.wilson_interval(5, 20, confidence=0.8)
    assert narrow[0] > 0.11186170140766569
    assert narrow[1] < 0.468700877618744


@pytest.mark.parametrize("confidence", [0, 1, -0.1, np.nan])
def test_invalid_confidence(confidence):
    with pytest.raises(ValueError, match="confidence"):
        cli.wilson_interval(1, 10, confidence)


@pytest.fixture
def output_dir(tmp_path):
    np.savetxt(tmp_path / "confined_fraction.dat", [
        [0, 0.5, 0.5, 10], [1e-4, 0.4, 0.4, 10], [1, 0.25, 0.25, 8]])
    return tmp_path


def test_losses_resolved_counts_and_log_axis(output_dir):
    destination = output_dir / "losses.png"
    figure = cli.plot_losses(output_dir, destination)
    assert destination.stat().st_size > 0
    ax = figure.axes[0]
    assert ax.get_xscale() == "log"
    for line, expected in zip(ax.lines, ([0.8, 0.5], [0.4, 0.25],
                                         [0.4, 0.25], [0.2, 0.5])):
        np.testing.assert_allclose(line.get_ydata(), expected)
        np.testing.assert_allclose(line.get_xdata(), [1e-4, 1])
    # Polygon bounds check the actual plotted band, including changing n.
    vertices = ax.collections[0].get_paths()[0].vertices
    for t, k, n in ((1e-4, 2, 10), (1, 4, 8)):
        bounds = vertices[vertices[:, 0] == t, 1]
        low, high = cli.wilson_interval(k, n)
        assert bounds.min() == pytest.approx(low)
        assert bounds.max() == pytest.approx(high)


def test_losses_legacy_marker_count_and_unresolved_gap(tmp_path):
    np.savetxt(tmp_path / "confined_fraction.dat", [[1e-4, 0.4, 0.4]])
    np.savetxt(tmp_path / "times_lost.dat", np.column_stack((np.arange(10), np.ones(10))))
    fig = cli.plot_losses(tmp_path, show=False)
    assert fig.axes[0].lines[0].get_ydata()[0] == pytest.approx(0.8)
    np.savetxt(tmp_path / "confined_fraction.dat", [
        [1e-4, 0.4, 0.4, 10], [1e-3, np.nan, np.nan, 0], [1, 0.25, 0.25, 8]])
    fig = cli.plot_losses(tmp_path, show=False)
    assert np.isnan(fig.axes[0].lines[0].get_ydata()[1])
    assert np.isnan(fig.axes[0].lines[3].get_ydata()[1])


@pytest.mark.parametrize("row", [[np.nan, 0.4, 0.4, 10], [1, 0.8, 0.8, 10], [1, -0.1, 0.5, 10],
                                   [1, 0.4, 0.4, 2.5], [1, np.nan, 0.4, 10]])
def test_losses_rejects_invalid_statistics(tmp_path, row):
    np.savetxt(tmp_path / "confined_fraction.dat", [row])
    with pytest.raises(ValueError, match="confined_fraction.dat"):
        cli.plot_losses(tmp_path, show=False)


def test_classes_text_counts(tmp_path):
    np.savetxt(tmp_path / "class_parts.dat", [
        [1, 0.2, 1, 1, 2, 0], [2, 0.2, 2, 2, 2, 1], [3, 0.3, 3, 1, 1, 2]])
    fig = cli.plot_classes(tmp_path, tmp_path / "classes.png")
    for ax, expected in zip(fig.axes, ([2, 1], [1, 2], [1, 1, 1])):
        assert [bar.get_height() for bar in ax.patches] == expected


@pytest.fixture(params=[False, True], ids=["native-dimensions", "transposed-dimensions"])
def results(tmp_path, request):
    nc = pytest.importorskip("netCDF4")
    path = tmp_path / "results.nc"
    with nc.Dataset(path, "w") as ds:
        for dim, n in (("particle", 3), ("xyz", 3), ("phase", 5), ("iclass_dim", 3)):
            ds.createDimension(dim, n)
        ds.cart_units = "cm"
        def write(name, data, component=None):
            dims = ("particle", component) if component else ("particle",)
            if component and request.param:
                dims = dims[::-1]
                data = np.asarray(data).T
            ds.createVariable(name, "f8", dims)[:] = data
        write("wall_hit", [1, 0, 1])
        write("wall_hit_cart", [[3, 4, 2], [900, 900, 900], [0, 2, -1]], "xyz")
        write("zend", [[0, 0, 0, 1, 0], [0, 0, 0, 99, 0], [0, 0, 0, 0.5, 0]], "phase")
        write("wall_hit_cos_incidence", [0.2, 0, 0.8])
        write("iclass", [[1, 2, 0], [2, 2, 1], [1, 1, 2]], "iclass_dim")
    return path


@pytest.mark.parametrize("color, expected", [("incidence", [0.2, 0.8]), ("energy", [1, 0.25])])
def test_wall_coordinates_mask_and_colour(results, color, expected):
    fig = cli.plot_wall(results, results.parent / f"wall-{color}.png", color_by=color)
    ax = fig.axes[0]
    np.testing.assert_allclose(ax.collections[0].get_offsets(), [[5, 2], [2, -1]])
    np.testing.assert_allclose(ax.collections[0].get_array(), expected)
    assert ax.get_xlabel() == "R / cm"
    assert ax.get_aspect() == 1


def test_classes_netcdf_named_dimensions(results):
    fig = cli.plot_classes(results, show=False)
    assert [[p.get_height() for p in ax.patches] for ax in fig.axes] == [
        [2, 1], [1, 2], [1, 1, 1]]


def test_wall_requires_wall_hits_not_generic_losses(results):
    from netCDF4 import Dataset
    with Dataset(results, "a") as ds:
        ds.renameVariable("wall_hit", "class_lost")
    with pytest.raises(ValueError, match="no wall_hit data"):
        cli.plot_wall(results, show=False)


def test_wall_no_hits(results):
    from netCDF4 import Dataset
    with Dataset(results, "a") as ds:
        ds.variables["wall_hit"][:] = 0
    with pytest.raises(ValueError, match="no wall-hit markers"):
        cli.plot_wall(results, show=False)


def test_poincare_coordinates_and_shared_colour_scale(tmp_path):
    for volume, s in ((1, -0.5), (2, 0.75)):
        np.savetxt(tmp_path / f"poincare_vol{volume}.dat",
                   [[1, 0, s, 0.1, 100, 10], [1, 1, s, 0.2, 101, 11]],
                   header="seed section s theta R Z")
    fig = cli.plot_poincare(tmp_path, tmp_path / "poincare.png")
    ax = fig.axes[0]
    np.testing.assert_allclose(ax.collections[0].get_offsets(), [[100, 10], [101, 11]])
    assert ax.collections[0].norm is ax.collections[1].norm
    assert (ax.collections[0].norm.vmin, ax.collections[0].norm.vmax) == (-0.5, 0.75)


@pytest.mark.parametrize("command, expected", [("losses", "confined_fraction.dat"),
    ("classes", "results.nc"), ("wall", "results.nc"), ("poincare", "poincare_vol")])
def test_cli_missing_output_is_clear(command, expected, tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main([command, str(tmp_path), "-o", str(tmp_path / "plot.png")])
    assert exc.value.code == 2
    error = capsys.readouterr().err
    assert expected in error
    assert "Traceback" not in error


def test_cli_dispatch_saves_without_display(output_dir, monkeypatch):
    def no_show():
        pytest.fail("saving a figure must not open an interactive window")
    monkeypatch.setattr(plt, "show", no_show)
    path = output_dir / "plot.png"
    assert cli.main(["losses", str(output_dir), "-o", str(path)]) == 0
    assert path.stat().st_size > 0
