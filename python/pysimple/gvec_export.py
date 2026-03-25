from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

import gvec

TESLA_TO_GAUSS = 1.0e4
METER_TO_CM = 1.0e2
WEBER_TO_GAUSS_CM2 = TESLA_TO_GAUSS * METER_TO_CM**2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a GVEC equilibrium to SIMPLE's NetCDF interchange format."
    )
    parser.add_argument("--param-file", required=True, type=Path)
    parser.add_argument("--state-file", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--family", choices=["logical"], default="logical")
    parser.add_argument("--ns", type=int, default=63)
    parser.add_argument("--ntheta", type=int, default=65)
    parser.add_argument("--nvarphi", type=int, default=65)
    parser.add_argument("--s-min", type=float, default=1.0e-12)
    parser.add_argument("--s-max", type=float, default=1.0)
    return parser.parse_args()


def _axis_scale(rho: np.ndarray) -> np.ndarray:
    scale = np.empty_like(rho)
    scale[0] = 1.0 / (2.0 * rho[0])
    scale[1:] = 1.0 / (2.0 * rho[1:])
    return scale


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


def _to_scalar_field(dataset: xr.Dataset, name: str) -> np.ndarray:
    return np.asarray(dataset[name].transpose("rad", "pol", "tor"), dtype=np.float64)


def _to_vector_field(dataset: xr.Dataset, name: str) -> np.ndarray:
    return np.asarray(
        dataset[name].transpose("rad", "pol", "tor", "xyz"), dtype=np.float64
    )


def _fortran_3d(values: np.ndarray) -> np.ndarray:
    return np.transpose(values, (2, 1, 0)).copy()


def build_dataset(args: argparse.Namespace) -> xr.Dataset:
    state = gvec.State(str(args.param_file), str(args.state_file))
    nfp = int(state.nfp)

    s = np.linspace(args.s_min, args.s_max, args.ns, dtype=np.float64)
    rho = np.sqrt(s)
    theta = np.linspace(0.0, 2.0 * np.pi, args.ntheta, dtype=np.float64)
    varphi = np.linspace(0.0, 2.0 * np.pi / nfp, args.nvarphi, dtype=np.float64)
    zeta = -varphi

    data = gvec.evaluate(
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

    phi = np.asarray(data["Phi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2
    chi = np.asarray(data["chi"], dtype=np.float64) * WEBER_TO_GAUSS_CM2
    iota = -np.asarray(data["iota"], dtype=np.float64)
    dphi_dr = np.asarray(data["dPhi_dr"], dtype=np.float64) * WEBER_TO_GAUSS_CM2

    scale = _axis_scale(rho)
    dphi_ds = dphi_dr * scale

    r_cm = _to_scalar_field(data, "X1") * METER_TO_CM
    z_cm = _to_scalar_field(data, "X2") * METER_TO_CM
    dr_ds_cm = _to_scalar_field(data, "dX1_dr") * scale[:, None, None] * METER_TO_CM
    dr_dt_cm = _to_scalar_field(data, "dX1_dt") * METER_TO_CM
    dr_dp_cm = -_to_scalar_field(data, "dX1_dz") * METER_TO_CM
    dz_ds_cm = _to_scalar_field(data, "dX2_dr") * scale[:, None, None] * METER_TO_CM
    dz_dt_cm = _to_scalar_field(data, "dX2_dt") * METER_TO_CM
    dz_dp_cm = -_to_scalar_field(data, "dX2_dz") * METER_TO_CM

    alam = _to_scalar_field(data, "LA")
    dl_ds = _to_scalar_field(data, "dLA_dr") * scale[:, None, None]
    dl_dt = _to_scalar_field(data, "dLA_dt")
    dl_dp = -_to_scalar_field(data, "dLA_dz")

    g_vmec, sqg_vmec = _metric_tensor_vmec(
        r_cm, dr_ds_cm, dr_dt_cm, dr_dp_cm, dz_ds_cm, dz_dt_cm, dz_dp_cm
    )
    g_sym, sqg = _metric_tensor_symflux(g_vmec, -sqg_vmec, dl_ds, dl_dt, dl_dp)

    dA_theta_ds = dphi_ds
    dA_phi_ds = -iota * dphi_ds
    bctr_vartheta = -dA_phi_ds[:, None, None] / sqg
    bctr_varphi = dA_theta_ds[:, None, None] / sqg

    bcov_s = g_sym[..., 0, 1] * bctr_vartheta + g_sym[..., 0, 2] * bctr_varphi
    bcov_vartheta = g_sym[..., 1, 1] * bctr_vartheta + g_sym[..., 1, 2] * bctr_varphi
    bcov_varphi = g_sym[..., 2, 1] * bctr_vartheta + g_sym[..., 2, 2] * bctr_varphi

    return xr.Dataset(
        data_vars={
            "A_theta": ("s", phi),
            "A_phi": ("s", chi),
            "dA_theta_ds": ("s", dA_theta_ds),
            "dA_phi_ds": ("s", dA_phi_ds),
            "iota": ("s", iota),
            "R": (("varphi", "theta", "s"), _fortran_3d(r_cm)),
            "Z": (("varphi", "theta", "s"), _fortran_3d(z_cm)),
            "dR_ds": (("varphi", "theta", "s"), _fortran_3d(dr_ds_cm)),
            "dR_dt": (("varphi", "theta", "s"), _fortran_3d(dr_dt_cm)),
            "dR_dp": (("varphi", "theta", "s"), _fortran_3d(dr_dp_cm)),
            "dZ_ds": (("varphi", "theta", "s"), _fortran_3d(dz_ds_cm)),
            "dZ_dt": (("varphi", "theta", "s"), _fortran_3d(dz_dt_cm)),
            "dZ_dp": (("varphi", "theta", "s"), _fortran_3d(dz_dp_cm)),
            "Lambda": (("varphi", "theta", "s"), _fortran_3d(alam)),
            "dLambda_ds": (("varphi", "theta", "s"), _fortran_3d(dl_ds)),
            "dLambda_dt": (("varphi", "theta", "s"), _fortran_3d(dl_dt)),
            "dLambda_dp": (("varphi", "theta", "s"), _fortran_3d(dl_dp)),
            "sqg_symflux": (("varphi", "theta", "s"), _fortran_3d(sqg)),
            "Bctr_vartheta": (("varphi", "theta", "s"), _fortran_3d(bctr_vartheta)),
            "Bctr_varphi": (("varphi", "theta", "s"), _fortran_3d(bctr_varphi)),
            "Bcov_s": (("varphi", "theta", "s"), _fortran_3d(bcov_s)),
            "Bcov_vartheta": (("varphi", "theta", "s"), _fortran_3d(bcov_vartheta)),
            "Bcov_varphi": (("varphi", "theta", "s"), _fortran_3d(bcov_varphi)),
        },
        coords={
            "s": s,
            "rho": ("s", rho),
            "theta": theta,
            "varphi": varphi,
        },
        attrs={
            "simple_file_type": "gvec_export",
            "simple_gvec_export_version": "1",
            "coordinate_family": args.family,
            "radial_coordinate": "s",
            "toroidal_coordinate": "varphi",
            "lambda_phi_convention": "simple_varphi",
            "length_unit": "cm",
            "vector_potential_unit": "gauss_cm2",
            "field_covariant_unit": "gauss_cm",
            "field_contravariant_unit": "1_per_cm",
            "jacobian_unit": "cm3",
            "nfp": nfp,
            "source_param_file": str(args.param_file),
            "source_state_file": str(args.state_file),
            "gvec_version": getattr(gvec, "__version__", "unknown"),
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
