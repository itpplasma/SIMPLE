// Direct-CUDA symplectic Euler (EXPL_IMPL_EULER) tracing on the Boozer
// splines, one particle per thread. The Newton driver mirrors gpu_newton1 and
// gpu_timestep_euler in simple_gpu.f90; every pointwise formula (spline cell
// polynomials, the guiding-centre chain, residual, Jacobian and the
// extrapolated angle advance) comes from gen/generate_orbit_systems.f90.
#include "simple_cuda_euler.h"

#include <cuda_runtime.h>
#include <math.h>
#include <stdint.h>

#include <chrono>

#include "generated/cuda_euler1_boozer.cuh"

namespace {

// FP64 throughput bounds the step on consumer GPUs: capping registers to raise
// occupancy (168 or 128 per thread) and blocks of 64 or 256 threads measured
// no faster than 128 threads at the compiler's 254 registers.
constexpr int kBlock = 128;

// Truncated constants as in libneo boozer_converter and SIMPLE util.
constexpr double kTwoPi = 6.28318530717958;
constexpr double kPi = 3.14159265358979;

enum {
  kStepOk = 0,
  kStepOutsideDomain = 1,
  kStepMaxIter = 2,
  kStepLinearSolve = 3,
  kStepBoundaryLimited = 6,
};
enum { kPredictNone = 0, kPredictLinear = 1 };
enum { kStopCorrection = 0 };
constexpr int kMaxIt = 32;

struct Spline1D {
  const double *coeff;
  int num_points;
  int periodic;
  double x_min, h_step;
};

struct Spline3D {
  const double *coeff;
  int num_points[3];
  int periodic[3];
  double x_min[3], h_step[3], inv_h_step[3], period[3];
};

struct Params {
  Spline1D aphi, bcov;
  Spline3D field;
  double torflux, phi_period;
  double rmax, stop_kappa;
  int predictor, stop_rule, warning_mode;
};

__device__ __forceinline__ double fmodulo(double a, double p) {
  return a - floor(a / p) * p;
}

// Cell location of libneo locate_cell_1d / locate_cell_3d.
__device__ __forceinline__ const double *cell_1d(const Spline1D &s, int nq,
                                                 double x, double *x_local) {
  double xj = x;
  if (s.periodic) {
    const double period = s.h_step * (double)(s.num_points - 1);
    xj = fmodulo(x - s.x_min, period) + s.x_min;
  }
  const double x_norm = (xj - s.x_min) / s.h_step;
  const int idx = max(0, min(s.num_points - 2, (int)x_norm));
  *x_local = (x_norm - (double)idx) * s.h_step;
  return s.coeff + (size_t)idx * (size_t)(6 * nq);
}

__device__ __forceinline__ const double *cell_3d(const Spline3D &s,
                                                 const double x[3],
                                                 double xl[3]) {
  int interval[3];
#pragma unroll
  for (int j = 0; j < 3; ++j) {
    double xj = x[j];
    if (s.periodic[j] &&
        (xj < s.x_min[j] || xj >= s.x_min[j] + s.period[j]))
      xj = fmodulo(xj - s.x_min[j], s.period[j]) + s.x_min[j];
    const double x_norm = (xj - s.x_min[j]) * s.inv_h_step[j];
    interval[j] = max(0, min(s.num_points[j] - 2, (int)x_norm));
    xl[j] = (x_norm - (double)interval[j]) * s.h_step[j];
  }
  const size_t cell =
      ((size_t)interval[2] * s.num_points[1] + interval[1]) * s.num_points[0] +
      interval[0];
  return s.coeff + cell * 216;
}

struct Aux {
  double pth, dpth_r, dpth_p, dh_r, d2h_rr, d2h_rp, d2pth_rr, d2pth_rp;
  double vpar, dvpar_r, dvpar_p, hth, dhth_r, hph, dhph_r;
};

// Field at (r, th, ph) with the axis continuation of gpu_eval_field_rmix
// (r < 0 maps to (|r|, th + pi)), then the generated Euler system.
__device__ __forceinline__ void euler_system(
    const Params &P, double r, double th, double ph, double pph, double mu,
    double ro0, double dt, double pthold, double z4, double fvec[2],
    double jac[4], Aux &a) {
  const bool crossed = r < 0.0;
  const double sigma = crossed ? -1.0 : 1.0;
  const double s = crossed ? -r : r;
  const double th_cont = crossed ? th + 0.5 * kTwoPi : th;

  double xl, aphi, aphi_r, aphi_rr;
  const double *c = cell_1d(P.aphi, 1, s, &xl);
  spline1d_o5_d2_nq1(c, xl, &aphi, &aphi_r, &aphi_rr);

  const double rho = sqrt(s);
  double x[3] = {rho, fmodulo(th_cont, kTwoPi), fmodulo(ph, P.phi_period)};
  double x3l[3], b, db[3], d2b[3];
  c = cell_3d(P.field, x, x3l);
  spline3d_o555_rmix_nq1(c, x3l[0], x3l[1], x3l[2], &b, db, d2b);

  double bc[2], dbc[2], d2bc[2];
  c = cell_1d(P.bcov, 2, rho, &xl);
  spline1d_o5_d2_nq2(c, xl, bc, dbc, d2bc);

  euler1_boozer_system(r, th, ph, pph, sigma, mu, ro0, P.torflux, dt, pthold,
                       z4, aphi, aphi_r, aphi_rr, bc[0], dbc[0], d2bc[0], bc[1],
                       dbc[1], d2bc[1], b, db[0], db[1], db[2], d2b[0], d2b[1],
                       d2b[2], &fvec[0], &fvec[1], &jac[0], &jac[1], &jac[2],
                       &jac[3], &a.pth, &a.dpth_r, &a.dpth_p, &a.dh_r,
                       &a.d2h_rr, &a.d2h_rp, &a.d2pth_rr, &a.d2pth_rp, &a.vpar,
                       &a.dvpar_r, &a.dvpar_p, &a.hth, &a.dhth_r, &a.hph,
                       &a.dhph_r);
}

struct State {
  double z[4], pthold, pth, dt, atol, rtol, warning_factor, mu, ro0;
  double hist_dt, hist_z[4], hist_dz[4][2];
  int ntau, nhist;
};

__device__ __forceinline__ bool history_valid(const State &S) {
  if (S.nhist <= 0) return false;
  return S.hist_dt == S.dt && S.hist_z[0] == S.z[0] && S.hist_z[1] == S.z[1] &&
         S.hist_z[2] == S.z[2] && S.hist_z[3] == S.z[3];
}

__device__ __forceinline__ void predict_state(const Params &P, const State &S,
                                              double zp[4]) {
  const double c = 1.0;
  for (int k = 0; k < 4; ++k) zp[k] = S.z[k];
  if (P.predictor == kPredictNone || !history_valid(S)) return;
  for (int k = 0; k < 4; ++k) zp[k] = S.z[k] + c * S.hist_dz[k][0];
  if (P.predictor != kPredictLinear && S.nhist > 1)
    for (int k = 0; k < 4; ++k)
      zp[k] = zp[k] + 0.5 * c * (c + 1.0) * (S.hist_dz[k][0] - S.hist_dz[k][1]);
  if (!(zp[0] > 0.0 && zp[0] < P.rmax))
    for (int k = 0; k < 4; ++k) zp[k] = S.z[k];
}

__device__ __forceinline__ void record_history(State &S, const double z0[4],
                                               bool keep) {
  if (!keep) {
    S.nhist = 0;
    return;
  }
  if (S.nhist > 0 && S.hist_dt == S.dt && S.hist_z[0] == z0[0] &&
      S.hist_z[1] == z0[1] && S.hist_z[2] == z0[2] && S.hist_z[3] == z0[3]) {
    for (int k = 0; k < 4; ++k) S.hist_dz[k][1] = S.hist_dz[k][0];
    S.nhist = min(S.nhist + 1, 2);
  } else {
    S.nhist = 1;
  }
  for (int k = 0; k < 4; ++k) {
    S.hist_dz[k][0] = S.z[k] - z0[k];
    S.hist_z[k] = S.z[k];
  }
  S.hist_dt = S.dt;
}

__device__ __forceinline__ bool rate_converged(const Params &P, double dnorm,
                                               double dnorm_prev) {
  if (P.stop_rule == kStopCorrection) return false;
  if (!(dnorm_prev > 0.0)) return false;
  const double theta = dnorm / dnorm_prev;
  if (!(theta < 1.0)) return false;
  return theta / (1.0 - theta) * dnorm <= P.stop_kappa;
}

__device__ __forceinline__ bool finite(double x) { return isfinite(x); }

// gpu_newton1 with sympl_euler1_newton_iter inlined.
__device__ int newton1(const Params &P, const State &S, double x[2],
                       double xlast[2], Aux &a, int *nfev) {
  int status = kStepMaxIter;
  *nfev = 0;
  bool boundary_limited = false, step_boundary_limited = false;
  double dnorm_prev = 0.0;
  const double tolref1 = 1.0;
  double tolref2 = fabs(1e1 * P.torflux / S.ro0);

  for (int kit = 1; kit <= kMaxIt; ++kit) {
    if (fabs(x[0]) > 1.0) return kStepOutsideDomain;
    double fvec[2], fjac[4];  // fjac column-major: 11, 21, 12, 22
    euler_system(P, x[0], S.z[1], S.z[2], x[1], S.mu, S.ro0, S.dt, S.pthold,
                 S.z[3], fvec, fjac, a);
    *nfev += 1;

    const double j11 = fjac[0], j21 = fjac[1], j12 = fjac[2], j22 = fjac[3];
    const double det = j11 * j22 - j12 * j21;
    const double scale =
        fmax(fmax(fabs(j11), fabs(j21)), fmax(fabs(j12), fabs(j22)));
    if (scale == 0.0 || fabs(det) <= 2.220446049250313e-16 * scale * scale) {
      xlast[0] = x[0];
      xlast[1] = x[1];
      return kStepLinearSolve;
    }
    double i11, i12, i21, i22;
    if (j11 != 0.0 && j21 != 0.0 && j12 != 0.0 && j22 != 0.0) {
      i11 = 1.0 / (j11 - j12 * j21 / j22);
      i12 = -1.0 / (j11 * j22 / j12 - j21);
      i21 = -1.0 / (j11 * j22 / j21 - j12);
      i22 = 1.0 / (j22 - j12 * j21 / j11);
    } else {
      i11 = j22 / det;
      i12 = -j12 / det;
      i21 = -j21 / det;
      i22 = j11 / det;
    }
    xlast[0] = x[0];
    xlast[1] = x[1];
    const double c1 = i11 * fvec[0] + i12 * fvec[1];
    const double c2 = i21 * fvec[0] + i22 * fvec[1];
    double step_scale = 1.0;
    if (P.rmax <= 1.0 && c1 < 0.0)
      step_scale = fmin(1.0, 0.8 * fmax(0.0, P.rmax - x[0]) / (-c1));
    step_boundary_limited = step_scale < 1.0;
    x[0] = x[0] - step_scale * c1;
    x[1] = x[1] - step_scale * c2;
    tolref2 = fmax(fabs(x[1]), tolref2);
    const bool converged =
        !step_boundary_limited &&
        ((fabs(fvec[0]) < S.atol && fabs(fvec[1]) < S.atol) ||
         (fabs(x[0] - xlast[0]) < S.rtol * tolref1 &&
          fabs(x[1] - xlast[1]) < S.rtol * tolref2));
    boundary_limited = boundary_limited || step_boundary_limited;

    if (fabs(x[0]) > 1.0) return kStepOutsideDomain;
    if (converged) return kStepOk;
    if (step_boundary_limited || !(S.rtol > 0.0)) {
      dnorm_prev = 0.0;
    } else {
      const double dnorm =
          fmax(fabs(x[0] - xlast[0]) / (S.rtol * tolref1),
               fabs(x[1] - xlast[1]) / (S.rtol * tolref2));
      if (rate_converged(P, dnorm, dnorm_prev)) return kStepOk;
      dnorm_prev = dnorm;
    }
  }
  if (boundary_limited) {
    if (step_boundary_limited &&
        fabs(1.0 - x[0]) <= fmax(1e-10, 10.0 * S.rtol))
      status = kStepOutsideDomain;
    else
      status = kStepBoundaryLimited;
  } else if (P.warning_mode && finite(x[0]) && finite(x[1]) &&
             finite(xlast[0]) && finite(xlast[1]) && finite(tolref2) &&
             fabs(x[0] - xlast[0]) <= S.warning_factor * S.rtol * tolref1 &&
             fabs(x[1] - xlast[1]) <= S.warning_factor * S.rtol * fabs(tolref2)) {
    status = kStepOk;
  }
  return status;
}

// gpu_timestep_euler.
__device__ int timestep_euler(const Params &P, State &S, int *nfev) {
  *nfev = 0;
  for (int ktau = 0; ktau < S.ntau; ++ktau) {
    const double accepted_pthold = S.pthold;
    S.pthold = S.pth;
    double z_start[4] = {S.z[0], S.z[1], S.z[2], S.z[3]};
    const bool predicted = history_valid(S);
    double zg[4];
    predict_state(P, S, zg);
    double x[2] = {zg[0], zg[3]}, xlast[2];
    Aux a;
    int n;
    int status = newton1(P, S, x, xlast, a, &n);
    *nfev += n;
    if (status != kStepOk && predicted) {
      S.pthold = accepted_pthold;
      x[0] = S.z[0];
      x[1] = S.z[3];
      status = newton1(P, S, x, xlast, a, &n);
      *nfev += n;
    }
    if (status != kStepOk) {
      S.pthold = accepted_pthold;
      return status;
    }
    if (x[0] > 1.0) {
      S.pthold = accepted_pthold;
      return kStepOutsideDomain;
    }
    bool crossed = false;
    if (x[0] < 0.0) {
      x[0] = -x[0];
      crossed = true;
      if (x[0] > 1.0) {
        S.pthold = accepted_pthold;
        return kStepOutsideDomain;
      }
      S.z[1] = S.z[1] + kPi;
    }
    S.z[0] = x[0];
    S.z[3] = x[1];
    double dr = x[0] - xlast[0], dp = x[1] - xlast[1];
    if (crossed) {
      // xlast lives in the other chart: evaluate the field afresh.
      double fvec[2], fjac[4];
      euler_system(P, S.z[0], S.z[1], S.z[2], S.z[3], S.mu, S.ro0, S.dt,
                   S.pthold, S.z[3], fvec, fjac, a);
      *nfev += 1;
      dr = 0.0;
      dp = 0.0;
    }
    euler1_boozer_advance(dr, dp, S.dt, S.z[1], S.z[2], a.pth, a.dpth_r,
                          a.dpth_p, a.dh_r, a.d2h_rr, a.d2h_rp, a.d2pth_rr,
                          a.d2pth_rp, a.vpar, a.dvpar_r, a.dvpar_p, a.hth,
                          a.dhth_r, a.hph, a.dhph_r, &S.pth, &S.z[1], &S.z[2]);
    record_history(S, z_start, !crossed);
  }
  return 0;
}

__global__ void __launch_bounds__(kBlock)
    trace_euler_kernel(Params P, int n, double *__restrict__ dstate,
                       int *__restrict__ istate, int ntimstep,
                       const int *__restrict__ ntau_macro,
                       int *__restrict__ loss_step,
                       double *__restrict__ loss_time,
                       int64_t *__restrict__ nfev_out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  State S;
  const double *d = dstate;
#define LD(field) S.field = d[(size_t)(k++) * n + i]
  int k = 0;
  for (int j = 0; j < 4; ++j) LD(z[j]);
  LD(pthold); LD(pth); LD(dt); LD(atol); LD(rtol); LD(warning_factor);
  LD(mu); LD(ro0); LD(hist_dt);
  for (int j = 0; j < 4; ++j) LD(hist_z[j]);
  for (int j = 0; j < 4; ++j) LD(hist_dz[j][0]);
  for (int j = 0; j < 4; ++j) LD(hist_dz[j][1]);
#undef LD
  S.ntau = istate[i];
  S.nhist = istate[n + i];

  int ierr = 0, lstep = ntimstep;
  int64_t nfev = 0;
  double elapsed = 0.0;
  for (int it = 2; it <= ntimstep && ierr == 0; ++it) {
    for (int kt = 0; kt < ntau_macro[it - 1]; ++kt) {
      int nf;
      ierr = timestep_euler(P, S, &nf);
      nfev += nf;
      elapsed += S.dt;
      if (ierr != 0) break;
    }
    if (ierr != 0) lstep = it;
  }
  if (ierr == 0) {
    int64_t total = 0;
    for (int it = 2; it <= ntimstep; ++it) total += ntau_macro[it - 1];
    elapsed = (double)total * S.dt;
  }
  loss_step[i] = lstep;
  loss_time[i] = elapsed;
  nfev_out[i] = nfev;

  double *o = dstate;
#define ST(field) o[(size_t)(k++) * n + i] = S.field
  k = 0;
  for (int j = 0; j < 4; ++j) ST(z[j]);
  ST(pthold); ST(pth); ST(dt); ST(atol); ST(rtol); ST(warning_factor);
  ST(mu); ST(ro0); ST(hist_dt);
  for (int j = 0; j < 4; ++j) ST(hist_z[j]);
  for (int j = 0; j < 4; ++j) ST(hist_dz[j][0]);
  for (int j = 0; j < 4; ++j) ST(hist_dz[j][1]);
#undef ST
  istate[n + i] = S.nhist;
}

__global__ void euler_system_probe_kernel(Params P, int n,
                                          const double *__restrict__ in,
                                          double *__restrict__ out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const double *p = in + (size_t)SIMPLE_CUDA_EULER_PROBE_INPUTS * i;
  double fvec[2], jac[4];
  Aux a;
  euler_system(P, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], fvec,
               jac, a);
  const double v[SIMPLE_CUDA_EULER_PROBE_OUTPUTS] = {
      fvec[0], fvec[1], jac[0], jac[1], jac[2], jac[3], a.pth, a.dpth_r,
      a.dpth_p, a.dh_r, a.d2h_rr, a.d2h_rp, a.d2pth_rr, a.d2pth_rp, a.vpar,
      a.dvpar_r, a.dvpar_p, a.hth, a.dhth_r, a.hph, a.dhph_r};
  for (int k = 0; k < SIMPLE_CUDA_EULER_PROBE_OUTPUTS; ++k)
    out[(size_t)SIMPLE_CUDA_EULER_PROBE_OUTPUTS * i + k] = v[k];
}

Params make_params(const struct simple_cuda_euler_field *field,
                   const struct simple_cuda_euler_controls *ctl,
                   double *const d_coeff[3]) {
  Params P{};
  P.aphi = {d_coeff[0], field->aphi_num_points, field->aphi_periodic,
            field->aphi_x_min, field->aphi_h_step};
  P.bcov = {d_coeff[1], field->bcov_num_points, field->bcov_periodic,
            field->bcov_x_min, field->bcov_h_step};
  P.field.coeff = d_coeff[2];
  for (int j = 0; j < 3; ++j) {
    P.field.num_points[j] = field->field_num_points[j];
    P.field.periodic[j] = field->field_periodic[j];
    P.field.x_min[j] = field->field_x_min[j];
    P.field.h_step[j] = field->field_h_step[j];
    P.field.inv_h_step[j] = field->field_inv_h_step[j];
    P.field.period[j] = field->field_period[j];
  }
  P.torflux = field->torflux;
  P.phi_period = kTwoPi / (double)field->nper;
  if (ctl) {
    P.rmax = ctl->rmax;
    P.stop_kappa = ctl->stop_kappa;
    P.predictor = ctl->predictor;
    P.stop_rule = ctl->stop_rule;
    P.warning_mode = ctl->warning_mode;
  }
  return P;
}

cudaError_t upload_field(const struct simple_cuda_euler_field *field,
                         double *d_coeff[3]) {
  const double *src[3] = {field->aphi_coeff, field->bcov_coeff,
                          field->field_coeff};
  const size_t count[3] = {field->aphi_count, field->bcov_count,
                           field->field_count};
  for (int k = 0; k < 3; ++k) {
    cudaError_t e = cudaMalloc(&d_coeff[k], count[k] * sizeof(double));
    if (e != cudaSuccess) return e;
    e = cudaMemcpy(d_coeff[k], src[k], count[k] * sizeof(double),
                   cudaMemcpyHostToDevice);
    if (e != cudaSuccess) return e;
  }
  return cudaSuccess;
}

double ms_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

}  // namespace

extern "C" int simple_cuda_euler_trace(
    int n, const struct simple_cuda_euler_field *field,
    const struct simple_cuda_euler_controls *ctl, double *state, int *istate,
    int ntimstep, const int *ntau_macro, int *loss_step, double *loss_time,
    int64_t *nfev, double profile_ms[SIMPLE_CUDA_EULER_PROFILE_COUNT]) {
  const auto t_total = std::chrono::steady_clock::now();
  for (int k = 0; k < SIMPLE_CUDA_EULER_PROFILE_COUNT; ++k) profile_ms[k] = 0;
  if (n <= 0) return 0;
  Params P{};
  double *d_coeff[3] = {nullptr, nullptr, nullptr};
  double *d_state = nullptr;
  int *d_istate = nullptr, *d_ntau = nullptr, *d_loss = nullptr;
  double *d_time = nullptr;
  int64_t *d_nfev = nullptr;
  cudaError_t e = cudaSuccess;
  const size_t nd = (size_t)SIMPLE_CUDA_EULER_STATE_COUNT * n;
  auto t0 = std::chrono::steady_clock::now();
#define CK(call)                \
  do {                          \
    e = (call);                 \
    if (e != cudaSuccess) goto done; \
  } while (0)
  CK(cudaMalloc(&d_state, nd * sizeof(double)));
  CK(cudaMalloc(&d_istate, 2 * (size_t)n * sizeof(int)));
  CK(cudaMalloc(&d_ntau, ntimstep * sizeof(int)));
  CK(cudaMalloc(&d_loss, n * sizeof(int)));
  CK(cudaMalloc(&d_time, n * sizeof(double)));
  CK(cudaMalloc(&d_nfev, n * sizeof(int64_t)));
  profile_ms[SIMPLE_CUDA_EULER_PROFILE_ALLOCATE] = ms_since(t0);

  t0 = std::chrono::steady_clock::now();
  CK(upload_field(field, d_coeff));
  CK(cudaMemcpy(d_state, state, nd * sizeof(double), cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d_istate, istate, 2 * (size_t)n * sizeof(int),
                cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d_ntau, ntau_macro, ntimstep * sizeof(int),
                cudaMemcpyHostToDevice));
  profile_ms[SIMPLE_CUDA_EULER_PROFILE_UPLOAD] = ms_since(t0);

  P = make_params(field, ctl, d_coeff);

  {
    t0 = std::chrono::steady_clock::now();
    trace_euler_kernel<<<(n + kBlock - 1) / kBlock, kBlock>>>(
        P, n, d_state, d_istate, ntimstep, d_ntau, d_loss, d_time, d_nfev);
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());
    profile_ms[SIMPLE_CUDA_EULER_PROFILE_KERNEL] = ms_since(t0);
  }

  t0 = std::chrono::steady_clock::now();
  CK(cudaMemcpy(state, d_state, nd * sizeof(double), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(istate, d_istate, 2 * (size_t)n * sizeof(int),
                cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(loss_step, d_loss, n * sizeof(int), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(loss_time, d_time, n * sizeof(double), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(nfev, d_nfev, n * sizeof(int64_t), cudaMemcpyDeviceToHost));
  profile_ms[SIMPLE_CUDA_EULER_PROFILE_DOWNLOAD] = ms_since(t0);
#undef CK
done:
  for (double *p : d_coeff) cudaFree(p);
  cudaFree(d_state);
  cudaFree(d_istate);
  cudaFree(d_ntau);
  cudaFree(d_loss);
  cudaFree(d_time);
  cudaFree(d_nfev);
  profile_ms[SIMPLE_CUDA_EULER_PROFILE_TOTAL] = ms_since(t_total);
  return (int)e;
}

extern "C" int simple_cuda_euler_system_probe(
    int n, const struct simple_cuda_euler_field *field, const double *in,
    double *out) {
  if (n <= 0) return 0;
  double *d_coeff[3] = {nullptr, nullptr, nullptr};
  double *d_in = nullptr, *d_out = nullptr;
  const size_t nin = (size_t)SIMPLE_CUDA_EULER_PROBE_INPUTS * n;
  const size_t nout = (size_t)SIMPLE_CUDA_EULER_PROBE_OUTPUTS * n;
  cudaError_t e = upload_field(field, d_coeff);
  if (e == cudaSuccess) e = cudaMalloc(&d_in, nin * sizeof(double));
  if (e == cudaSuccess) e = cudaMalloc(&d_out, nout * sizeof(double));
  if (e == cudaSuccess)
    e = cudaMemcpy(d_in, in, nin * sizeof(double), cudaMemcpyHostToDevice);
  if (e == cudaSuccess) {
    euler_system_probe_kernel<<<(n + kBlock - 1) / kBlock, kBlock>>>(
        make_params(field, nullptr, d_coeff), n, d_in, d_out);
    e = cudaGetLastError();
  }
  if (e == cudaSuccess)
    e = cudaMemcpy(out, d_out, nout * sizeof(double), cudaMemcpyDeviceToHost);
  for (double *p : d_coeff) cudaFree(p);
  cudaFree(d_in);
  cudaFree(d_out);
  return (int)e;
}
