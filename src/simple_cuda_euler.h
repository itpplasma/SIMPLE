#ifndef SIMPLE_CUDA_EULER_H
#define SIMPLE_CUDA_EULER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  SIMPLE_CUDA_EULER_PROFILE_ALLOCATE = 0,
  SIMPLE_CUDA_EULER_PROFILE_UPLOAD = 1,
  SIMPLE_CUDA_EULER_PROFILE_KERNEL = 2,
  SIMPLE_CUDA_EULER_PROFILE_DOWNLOAD = 3,
  SIMPLE_CUDA_EULER_PROFILE_TOTAL = 4,
  SIMPLE_CUDA_EULER_PROFILE_COUNT = 5,
};

/* Per-particle real state, structure of arrays state[field*n + i] in this
   order: z(1:4), pthold, pth, dt, atol, rtol, warning_factor, mu, ro0,
   hist_dt, hist_z(1:4), hist_dz(1:4,1), hist_dz(1:4,2). The integer state is
   istate[i] = ntau, istate[n + i] = nhist. */
enum { SIMPLE_CUDA_EULER_STATE_COUNT = 25 };

/* Quintic libneo batch splines: A_phi(s) (one quantity), B_vartheta and
   B_varphi (rho) (two quantities) and |B|(rho, vartheta, varphi), with the
   coefficient arrays in libneo's Fortran layout. */
struct simple_cuda_euler_field {
  const double *aphi_coeff;
  size_t aphi_count;
  int aphi_num_points, aphi_periodic;
  double aphi_x_min, aphi_h_step;
  const double *bcov_coeff;
  size_t bcov_count;
  int bcov_num_points, bcov_periodic;
  double bcov_x_min, bcov_h_step;
  const double *field_coeff;
  size_t field_count;
  int field_num_points[3], field_periodic[3];
  double field_x_min[3], field_h_step[3], field_inv_h_step[3], field_period[3];
  double torflux;
  int nper;
};

struct simple_cuda_euler_controls {
  double rmax, stop_kappa;
  int predictor, stop_rule, warning_mode;
};

int simple_cuda_euler_trace(
    int n, const struct simple_cuda_euler_field *field,
    const struct simple_cuda_euler_controls *controls, double *state,
    int *istate, int ntimstep, const int *ntau_macro, int *loss_step,
    double *loss_time, int64_t *nfev,
    double profile_ms[SIMPLE_CUDA_EULER_PROFILE_COUNT]);

/* Field and generated Euler system at n points, for tests. Inputs per point:
   r, theta, phi, p_phi, mu, ro0, dt, pthold, z4. Outputs per point: fvec(2),
   the Jacobian (11, 21, 12, 22), p_theta, dp_theta/dr, dp_theta/dp_phi,
   dH/dr, d2H/dr2, d2H/drdp_phi, d2p_theta/dr2, d2p_theta/drdp_phi, v_par,
   dv_par/dr, dv_par/dp_phi, h_theta, dh_theta/dr, h_phi, dh_phi/dr. */
enum {
  SIMPLE_CUDA_EULER_PROBE_INPUTS = 9,
  SIMPLE_CUDA_EULER_PROBE_OUTPUTS = 21,
};

int simple_cuda_euler_system_probe(
    int n, const struct simple_cuda_euler_field *field, const double *in,
    double *out);

#ifdef __cplusplus
}
#endif

#endif
