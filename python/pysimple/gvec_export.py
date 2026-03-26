from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

import gvec

TESLA_TO_GAUSS = 1.0e4
METER_TO_CM = 1.0e2
WEBER_TO_GAUSS_CM2 = TESLA_TO_GAUSS * METER_TO_CM**2
SFL_RHO_MIN = 1.1e-4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a GVEC equilibrium to SIMPLE's NetCDF interchange format."
    )
    parser.add_argument("--param-file", required=True, type=Path)
    parser.add_argument("--state-file", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--family", choices=["logical", "boozer"], default="logical")
    parser.add_argument("--ns", type=int, default=63)
    parser.add_argument("--ntheta", type=int, default=65)
    parser.add_argument("--nvarphi", type=int, default=65)
    parser.add_argument("--s-min", type=float, default=1.0e-12)
    parser.add_argument("--s-max", type=float, default=1.0)
    parser.add_argument("--mnfactor", type=int, default=1)
    return parser.parse_args()


def _axis_scale(rho: np.ndarray) -> np.ndarray:
    scale = np.empty_like(rho)
    scale[0] = 1.0 / (2.0 * rho[0])
    scale[1:] = 1.0 / (2.0 * rho[1:])
    return scale


def _transpose_xyz(dataset: xr.Dataset, name: str) -> np.ndarray:
    return np.asarray(
        dataset[name].transpose("rad", "pol", "tor", "xyz"), dtype=np.float64
    )


def _transpose_scalar(dataset: xr.Dataset, name: str) -> np.ndarray:
    return np.asarray(dataset[name].transpose("rad", "pol", "tor"), dtype=np.float64)


def _fortran_1d(values: np.ndarray) -> np.ndarray:
    return np.asarray(values, dtype=np.float64).copy()


def _fortran_3d(values: np.ndarray) -> np.ndarray:
    return np.transpose(values, (2, 1, 0)).copy()


def _geometry_from_basis(
    pos_cm: np.ndarray,
    e_s_cm: np.ndarray,
    e_t_cm: np.ndarray,
    e_p_cm: np.ndarray,
) -> dict[str, np.ndarray]:
    return {
        "X": pos_cm[..., 0],
        "Y": pos_cm[..., 1],
        "Z": pos_cm[..., 2],
        "dX_ds": e_s_cm[..., 0],
        "dX_dt": e_t_cm[..., 0],
        "dX_dp": e_p_cm[..., 0],
        "dY_ds": e_s_cm[..., 1],
        "dY_dt": e_t_cm[..., 1],
        "dY_dp": e_p_cm[..., 1],
        "dZ_ds": e_s_cm[..., 2],
        "dZ_dt": e_t_cm[..., 2],
        "dZ_dp": e_p_cm[..., 2],
    }


def _broadcast_profile(profile: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    return np.broadcast_to(profile[:, None, None], shape).copy()


def _metric_tensor_vmec(
    r_cm: np.ndarray,
    dr_ds_cm: np.ndarray,
    dr_dt_cm: np.ndarray,
    dr_dp_cm: np.ndarray,
    dz_ds_cm: np.ndarray,
    dz_dt_cm: np.ndarray,
    dz_dp_cm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    g = np.zeros(r_cm.shape + (3, 3), dtype=np.float64)

    g[..., 0, 0] = dr_ds_cm**2 + dz_ds_cm**2
    g[..., 0, 1] = dr_ds_cm * dr_dt_cm + dz_ds_cm * dz_dt_cm
    g[..., 1, 0] = g[..., 0, 1]
    g[..., 0, 2] = dr_ds_cm * dr_dp_cm + dz_ds_cm * dz_dp_cm
    g[..., 2, 0] = g[..., 0, 2]
    g[..., 1, 1] = dr_dt_cm**2 + dz_dt_cm**2
    g[..., 1, 2] = dr_dt_cm * dr_dp_cm + dz_dt_cm * dz_dp_cm
    g[..., 2, 1] = g[..., 1, 2]
    g[..., 2, 2] = r_cm**2 + dr_dp_cm**2 + dz_dp_cm**2

    sqg_vmec = -r_cm * (dz_ds_cm * dr_dt_cm - dr_ds_cm * dz_dt_cm)
    return g, sqg_vmec


def _metric_tensor_symflux(
    g_vmec: np.ndarray,
    sqg_vmec: np.ndarray,
    dl_ds: np.ndarray,
    dl_dt: np.ndarray,
    dl_dp: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    cjac = 1.0 / (1.0 + dl_dt)
    cmat = np.zeros(g_vmec.shape, dtype=np.float64)
    cmat[..., 0, 0] = 1.0
    cmat[..., 1, 0] = -dl_ds * cjac
    cmat[..., 1, 1] = cjac
    cmat[..., 1, 2] = -dl_dp * cjac
    cmat[..., 2, 2] = 1.0
    g_sym = np.einsum("...ji,...jk,...kl->...il", cmat, g_vmec, cmat)
    sqg = sqg_vmec * cjac
    return g_sym, sqg


def _jacobian_from_basis(
    e_s_cm: np.ndarray,
    e_t_cm: np.ndarray,
    e_p_cm: np.ndarray,
) -> np.ndarray:
    return np.einsum("...i,...i->...", e_s_cm, np.cross(e_t_cm, e_p_cm, axis=-1))


def _covariant_from_vector(
    vector_xyz: np.ndarray,
    e_s_cm: np.ndarray,
    e_t_cm: np.ndarray,
    e_p_cm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cov_s = np.einsum("...i,...i->...", vector_xyz, e_s_cm)
    cov_t = np.einsum("...i,...i->...", vector_xyz, e_t_cm)
    cov_p = np.einsum("...i,...i->...", vector_xyz, e_p_cm)
    return cov_s, cov_t, cov_p


def _logical_dataset(
    state: gvec.State,
    rho: np.ndarray,
    theta: np.ndarray,
    varphi: np.ndarray,
) -> xr.Dataset:
    zeta = -varphi
    return gvec.evaluate(
        state,
        "Phi",
        "chi",
        "iota",
        "LA",
        "dLA_dr",
        "dLA_dt",
        "dLA_dz",
        "dPhi_dr",
        "X1",
        "X2",
        "dX1_dr",
        "dX1_dt",
        "dX1_dz",
        "dX2_dr",
        "dX2_dt",
        "dX2_dz",
        rho=rho,
        theta=theta,
        zeta=zeta,
    )


def _boozer_dataset(
    state: gvec.State,
    rho: np.ndarray,
    theta: np.ndarray,
    varphi: np.ndarray,
    mnfactor: int,
) -> xr.Dataset:
    return gvec.evaluate_sfl(
        state,
        "pos",
        "B",
        "mod_B",
        "e_rho_B",
        "e_theta_B",
        "e_zeta_B",
        "Phi",
        "dPhi_dr",
        "chi",
        "dchi_dr",
        "iota",
        "B_rho_B",
        "B_theta_B",
        "B_zeta_B",
        "B_contra_t_B",
        "B_contra_z_B",
        rho=rho,
        theta=theta,
        zeta=varphi,
        sfl="boozer",
        MNfactor=mnfactor,
    )


def _logical_export(
    ds: xr.Dataset,
    rho: np.ndarray,
    scale_1d: np.ndarray,
    scale_3d: np.ndarray,
    phi_period: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, str]]:
    phi = _fortran_1d(np.asarray(ds["Phi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2)
    chi = _fortran_1d(np.asarray(ds["chi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2)
    iota = _fortran_1d(-np.asarray(ds["iota"], dtype=np.float64))
    dphi_ds = _fortran_1d(
        np.asarray(ds["dPhi_dr"], dtype=np.float64) * scale_1d * WEBER_TO_GAUSS_CM2
    )
    dchi_ds = -iota * dphi_ds

    r_cm = _transpose_scalar(ds, "X1") * METER_TO_CM
    z_cm = _transpose_scalar(ds, "X2") * METER_TO_CM
    dr_ds_cm = _transpose_scalar(ds, "dX1_dr") * scale_3d[..., 0] * METER_TO_CM
    dr_dt_cm = _transpose_scalar(ds, "dX1_dt") * METER_TO_CM
    dr_dp_cm = -_transpose_scalar(ds, "dX1_dz") * METER_TO_CM
    dz_ds_cm = _transpose_scalar(ds, "dX2_dr") * scale_3d[..., 0] * METER_TO_CM
    dz_dt_cm = _transpose_scalar(ds, "dX2_dt") * METER_TO_CM
    dz_dp_cm = -_transpose_scalar(ds, "dX2_dz") * METER_TO_CM

    lambda_vals = _transpose_scalar(ds, "LA")
    dl_ds = _transpose_scalar(ds, "dLA_dr") * scale_3d[..., 0]
    dl_dt = _transpose_scalar(ds, "dLA_dt")
    dl_dp = -_transpose_scalar(ds, "dLA_dz")

    g_vmec, sqg_vmec = _metric_tensor_vmec(
        r_cm, dr_ds_cm, dr_dt_cm, dr_dp_cm, dz_ds_cm, dz_dt_cm, dz_dp_cm
    )
    g_sym, sqg = _metric_tensor_symflux(g_vmec, -sqg_vmec, dl_ds, dl_dt, dl_dp)

    bctr_t = -dchi_ds[:, None, None] / sqg
    bctr_p = dphi_ds[:, None, None] / sqg
    bcov_s = g_sym[..., 0, 1] * bctr_t + g_sym[..., 0, 2] * bctr_p
    bcov_t = g_sym[..., 1, 1] * bctr_t + g_sym[..., 1, 2] * bctr_p
    bcov_p = g_sym[..., 2, 1] * bctr_t + g_sym[..., 2, 2] * bctr_p
    bmod = np.sqrt(bctr_t * bcov_t + bctr_p * bcov_p)
    grid_shape = bmod.shape

    varphi = np.linspace(0.0, phi_period, grid_shape[2], endpoint=False, dtype=np.float64)
    cos_phi = np.cos(varphi)[None, None, :]
    sin_phi = np.sin(varphi)[None, None, :]
    geometry = {
        "R": r_cm,
        "Z_cyl": z_cm,
        "dR_ds": dr_ds_cm,
        "dR_dt": dr_dt_cm,
        "dR_dp": dr_dp_cm,
        "dZ_ds_cyl": dz_ds_cm,
        "dZ_dt_cyl": dz_dt_cm,
        "dZ_dp_cyl": dz_dp_cm,
        "X": r_cm * cos_phi,
        "Y": r_cm * sin_phi,
        "Z": z_cm,
        "dX_ds": dr_ds_cm * cos_phi,
        "dX_dt": dr_dt_cm * cos_phi,
        "dX_dp": dr_dp_cm * cos_phi - r_cm * sin_phi,
        "dY_ds": dr_ds_cm * sin_phi,
        "dY_dt": dr_dt_cm * sin_phi,
        "dY_dp": dr_dp_cm * sin_phi + r_cm * cos_phi,
        "dZ_ds": dz_ds_cm,
        "dZ_dt": dz_dt_cm,
        "dZ_dp": dz_dp_cm,
    }

    a_theta = _broadcast_profile(phi, grid_shape)
    a_phi = _broadcast_profile(chi, grid_shape)
    acov_s = a_theta * dl_ds
    acov_t = a_theta * (1.0 + dl_dt)
    acov_p = a_phi + a_theta * dl_dp
    hcov_s = bcov_s / bmod
    hcov_t = bcov_t / bmod
    hcov_p = bcov_p / bmod

    field = {
        "Acov_s": acov_s,
        "Acov_t": acov_t,
        "Acov_p": acov_p,
        "hcov_s": hcov_s,
        "hcov_t": hcov_t,
        "hcov_p": hcov_p,
        "Bmod": bmod,
        "sqgBctr_t": _broadcast_profile(-dchi_ds, grid_shape),
        "sqgBctr_p": _broadcast_profile(dphi_ds, grid_shape),
    }
    profiles = {
        "A_theta": phi,
        "A_phi": chi,
        "dA_theta_ds": dphi_ds,
        "dA_phi_ds": dchi_ds,
        "iota": iota,
    }
    transforms = {
        "Lambda": lambda_vals,
        "dLambda_ds": dl_ds,
        "dLambda_dt": dl_dt,
        "dLambda_dp": dl_dp,
    }
    attrs = {
        "coordinate_family": "logical",
        "theta_coordinate_name": "theta",
        "phi_coordinate_name": "varphi",
        "phi_period": str(phi_period),
        "has_lambda_transform": "true",
    }
    return geometry | field | profiles | transforms, profiles, attrs


def _straight_export(
    ds: xr.Dataset,
    rho: np.ndarray,
    scale_1d: np.ndarray,
    scale_3d: np.ndarray,
    family: str,
    phi_period: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, str]]:
    pos_cm = _transpose_xyz(ds, "pos") * METER_TO_CM
    b_xyz = _transpose_xyz(ds, "B") * TESLA_TO_GAUSS

    basis_prefix = "B" if family == "boozer" else "P"
    e_s_cm = _transpose_xyz(ds, f"e_rho_{basis_prefix}") * scale_3d * METER_TO_CM
    e_t_cm = _transpose_xyz(ds, f"e_theta_{basis_prefix}") * METER_TO_CM
    e_p_cm = _transpose_xyz(ds, f"e_zeta_{basis_prefix}") * METER_TO_CM

    phi = _fortran_1d(np.asarray(ds["Phi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2)
    chi = _fortran_1d(np.asarray(ds["chi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2)
    iota = _fortran_1d(np.asarray(ds["iota"], dtype=np.float64))
    dphi_ds = _fortran_1d(
        np.asarray(ds["dPhi_dr"], dtype=np.float64) * scale_1d * WEBER_TO_GAUSS_CM2
    )
    dchi_ds = _fortran_1d(
        np.asarray(ds["dchi_dr"], dtype=np.float64) * scale_1d * WEBER_TO_GAUSS_CM2
    )

    geometry = _geometry_from_basis(pos_cm, e_s_cm, e_t_cm, e_p_cm)
    sqg = _jacobian_from_basis(e_s_cm, e_t_cm, e_p_cm)
    bmod = np.linalg.norm(b_xyz, axis=-1)
    grid_shape = bmod.shape
    bcov_rho = (
        _transpose_scalar(ds, f"B_rho_{basis_prefix}") * TESLA_TO_GAUSS * METER_TO_CM
    )
    bcov_t = (
        _transpose_scalar(ds, f"B_theta_{basis_prefix}") * TESLA_TO_GAUSS * METER_TO_CM
    )
    bcov_p = (
        _transpose_scalar(ds, f"B_zeta_{basis_prefix}") * TESLA_TO_GAUSS * METER_TO_CM
    )
    bctr_t = _transpose_scalar(ds, f"B_contra_t_{basis_prefix}") * (
        TESLA_TO_GAUSS / METER_TO_CM
    )
    bctr_p = _transpose_scalar(ds, f"B_contra_z_{basis_prefix}") * (
        TESLA_TO_GAUSS / METER_TO_CM
    )
    bcov_s = bcov_rho * scale_3d[..., 0]
    sqgbctr_t = sqg * bctr_t
    sqgbctr_p = sqg * bctr_p

    acov_s = np.zeros_like(bmod)
    acov_t = _broadcast_profile(phi, grid_shape)
    acov_p = _broadcast_profile(chi, grid_shape)
    hcov_s = bcov_s / bmod
    hcov_t = bcov_t / bmod
    hcov_p = bcov_p / bmod

    field = {
        "Acov_s": acov_s,
        "Acov_t": acov_t,
        "Acov_p": acov_p,
        "hcov_s": hcov_s,
        "hcov_t": hcov_t,
        "hcov_p": hcov_p,
        "Bmod": bmod,
        "sqgBctr_t": sqgbctr_t,
        "sqgBctr_p": sqgbctr_p,
        "B_x": b_xyz[..., 0],
        "B_y": b_xyz[..., 1],
        "B_z": b_xyz[..., 2],
    }
    profiles = {
        "A_theta": phi,
        "A_phi": chi,
        "dA_theta_ds": dphi_ds,
        "dA_phi_ds": dchi_ds,
        "iota": iota,
    }
    attrs = {
        "coordinate_family": family,
        "theta_coordinate_name": "theta_B",
        "phi_coordinate_name": "zeta_B",
        "phi_period": str(phi_period),
        "has_lambda_transform": "false",
    }
    return geometry | field | profiles, profiles, attrs


def build_dataset(args: argparse.Namespace) -> xr.Dataset:
    state = gvec.State(str(args.param_file), str(args.state_file))
    nfp = int(state.nfp)

    if args.family == "logical":
        s_min = args.s_min
    else:
        s_min = max(args.s_min, SFL_RHO_MIN**2)

    s = np.linspace(s_min, args.s_max, args.ns, dtype=np.float64)
    rho = np.sqrt(s)
    theta_endpoint = args.family == "logical"
    varphi_endpoint = args.family == "logical"
    theta = np.linspace(
        0.0, 2.0 * np.pi, args.ntheta, endpoint=theta_endpoint, dtype=np.float64
    )
    varphi = np.linspace(
        0.0,
        2.0 * np.pi / nfp,
        args.nvarphi,
        endpoint=varphi_endpoint,
        dtype=np.float64,
    )
    scale_1d = _axis_scale(rho)
    scale_3d = scale_1d[:, None, None, None]
    phi_period = 2.0 * np.pi / nfp

    if args.family == "logical":
        raw = _logical_dataset(state, rho, theta, varphi)
        data, _, family_attrs = _logical_export(
            raw, rho, scale_1d, scale_3d, phi_period
        )
    elif args.family == "boozer":
        raw = _boozer_dataset(state, rho, theta, varphi, args.mnfactor)
        data, _, family_attrs = _straight_export(
            raw, rho, scale_1d, scale_3d, args.family, phi_period
        )
    else:
        raise ValueError(f"Unsupported export family: {args.family}")

    data_vars = {
        "A_theta": ("s", _fortran_1d(data["A_theta"])),
        "A_phi": ("s", _fortran_1d(data["A_phi"])),
        "dA_theta_ds": ("s", _fortran_1d(data["dA_theta_ds"])),
        "dA_phi_ds": ("s", _fortran_1d(data["dA_phi_ds"])),
        "iota": ("s", _fortran_1d(data["iota"])),
        "X": (("varphi", "theta", "s"), _fortran_3d(data["X"])),
        "Y": (("varphi", "theta", "s"), _fortran_3d(data["Y"])),
        "Z": (("varphi", "theta", "s"), _fortran_3d(data["Z"])),
        "dX_ds": (("varphi", "theta", "s"), _fortran_3d(data["dX_ds"])),
        "dX_dt": (("varphi", "theta", "s"), _fortran_3d(data["dX_dt"])),
        "dX_dp": (("varphi", "theta", "s"), _fortran_3d(data["dX_dp"])),
        "dY_ds": (("varphi", "theta", "s"), _fortran_3d(data["dY_ds"])),
        "dY_dt": (("varphi", "theta", "s"), _fortran_3d(data["dY_dt"])),
        "dY_dp": (("varphi", "theta", "s"), _fortran_3d(data["dY_dp"])),
        "dZ_ds": (("varphi", "theta", "s"), _fortran_3d(data["dZ_ds"])),
        "dZ_dt": (("varphi", "theta", "s"), _fortran_3d(data["dZ_dt"])),
        "dZ_dp": (("varphi", "theta", "s"), _fortran_3d(data["dZ_dp"])),
        "Acov_s": (("varphi", "theta", "s"), _fortran_3d(data["Acov_s"])),
        "Acov_t": (("varphi", "theta", "s"), _fortran_3d(data["Acov_t"])),
        "Acov_p": (("varphi", "theta", "s"), _fortran_3d(data["Acov_p"])),
        "hcov_s": (("varphi", "theta", "s"), _fortran_3d(data["hcov_s"])),
        "hcov_t": (("varphi", "theta", "s"), _fortran_3d(data["hcov_t"])),
        "hcov_p": (("varphi", "theta", "s"), _fortran_3d(data["hcov_p"])),
        "Bmod": (("varphi", "theta", "s"), _fortran_3d(data["Bmod"])),
        "sqgBctr_t": (("varphi", "theta", "s"), _fortran_3d(data["sqgBctr_t"])),
        "sqgBctr_p": (("varphi", "theta", "s"), _fortran_3d(data["sqgBctr_p"])),
    }
    if "R" in data:
        data_vars["R"] = (("varphi", "theta", "s"), _fortran_3d(data["R"]))
        data_vars["Z_cyl"] = (("varphi", "theta", "s"), _fortran_3d(data["Z_cyl"]))
        data_vars["dR_ds"] = (("varphi", "theta", "s"), _fortran_3d(data["dR_ds"]))
        data_vars["dR_dt"] = (("varphi", "theta", "s"), _fortran_3d(data["dR_dt"]))
        data_vars["dR_dp"] = (("varphi", "theta", "s"), _fortran_3d(data["dR_dp"]))
        data_vars["dZ_ds_cyl"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dZ_ds_cyl"]),
        )
        data_vars["dZ_dt_cyl"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dZ_dt_cyl"]),
        )
        data_vars["dZ_dp_cyl"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dZ_dp_cyl"]),
        )
    if "Lambda" in data:
        data_vars["Lambda"] = (("varphi", "theta", "s"), _fortran_3d(data["Lambda"]))
        data_vars["dLambda_ds"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dLambda_ds"]),
        )
        data_vars["dLambda_dt"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dLambda_dt"]),
        )
        data_vars["dLambda_dp"] = (
            ("varphi", "theta", "s"),
            _fortran_3d(data["dLambda_dp"]),
        )
    if "B_x" in data:
        data_vars["B_x"] = (("varphi", "theta", "s"), _fortran_3d(data["B_x"]))
        data_vars["B_y"] = (("varphi", "theta", "s"), _fortran_3d(data["B_y"]))
        data_vars["B_z"] = (("varphi", "theta", "s"), _fortran_3d(data["B_z"]))

    return xr.Dataset(
        data_vars=data_vars,
        coords={
            "s": s,
            "rho": ("s", rho),
            "theta": theta,
            "varphi": varphi,
        },
        attrs={
            "simple_file_type": "gvec_export",
            "simple_gvec_export_version": "1",
            "radial_coordinate": "s",
            "length_unit": "cm",
            "vector_potential_unit": "gauss_cm2",
            "field_covariant_unit": "gauss_cm",
            "field_direction_unit": "1",
            "sqgBctr_unit": "gauss_cm",
            "nfp": nfp,
            "source_param_file": str(args.param_file),
            "source_state_file": str(args.state_file),
            "gvec_version": getattr(gvec, "__version__", "unknown"),
            **family_attrs,
        },
    )


def main() -> None:
    args = parse_args()
    dataset = build_dataset(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = args.output.with_suffix(args.output.suffix + ".tmp")
    if tmp_output.exists():
        tmp_output.unlink()
    dataset.to_netcdf(tmp_output)
    tmp_output.replace(args.output)


if __name__ == "__main__":
    main()
