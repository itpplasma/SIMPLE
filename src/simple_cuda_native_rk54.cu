#include "simple_cuda_native_rk54.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>

namespace {

#ifndef SIMPLE_CUDA_NATIVE_THREADS
#define SIMPLE_CUDA_NATIVE_THREADS 32
#endif

constexpr int kThreads = SIMPLE_CUDA_NATIVE_THREADS;
constexpr int kMaximumAttempts = 100000000;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kSqrt2 = 1.414213562373095048801688724209698079;

static_assert(kThreads > 0 && kThreads <= 1024,
              "SIMPLE_CUDA_NATIVE_THREADS must be in [1, 1024]");

struct Geometry {
  int point_count[3];
  double x_min[3];
  double inv_h_step[3];
  double period[3];
  double inv_period[3];
  double torflux;
  double total_duration;
  double block_duration;
  double tolerance;
  double minimum_timestep;
};

__device__ __forceinline__ void
cubic_table_location(double x, int dimension, const Geometry &geometry,
                     int &first, float weight[4], float *derivative) {
  double x_eval = x;
  if (dimension > 0) {
    const double periods =
        floor((x - geometry.x_min[dimension]) * geometry.inv_period[dimension]);
    x_eval = x - periods * geometry.period[dimension];
  }
  const double x_grid =
      (x_eval - geometry.x_min[dimension]) * geometry.inv_h_step[dimension];
  int first_zero = 3 * (static_cast<int>(x_grid) / 3);
  first_zero = max(0, min(first_zero, geometry.point_count[dimension] - 4));
  const float relative = static_cast<float>(x_grid - first_zero);
  first = first_zero;

  weight[0] = (1.0f - relative) * (2.0f - relative) * (3.0f - relative) / 6.0f;
  weight[1] = relative * (2.0f - relative) * (3.0f - relative) / 2.0f;
  weight[2] = relative * (relative - 1.0f) * (3.0f - relative) / 2.0f;
  weight[3] = relative * (relative - 1.0f) * (relative - 2.0f) / 6.0f;
  if (derivative) {
    derivative[0] = (-11.0f + relative * (12.0f - 3.0f * relative)) / 6.0f;
    derivative[1] = (6.0f + relative * (-10.0f + 3.0f * relative)) / 2.0f;
    derivative[2] = (-3.0f + relative * (8.0f - 3.0f * relative)) / 2.0f;
    derivative[3] = (2.0f + relative * (-6.0f + 3.0f * relative)) / 6.0f;
#pragma unroll
    for (int i = 0; i < 4; ++i)
      derivative[i] *= static_cast<float>(geometry.inv_h_step[dimension]);
  }
}

__device__ __forceinline__ void
evaluate_profile_table(const float *__restrict__ table,
                       const Geometry &geometry, double s, float value[6]) {
  int first;
  float weight[4];
  cubic_table_location(s, 0, geometry, first, weight, nullptr);
#pragma unroll
  for (int q = 0; q < 6; ++q)
    value[q] = 0.0f;
#pragma unroll
  for (int i = 0; i < 4; ++i) {
    const float *point = table + 6 * (first + i);
#pragma unroll
    for (int q = 0; q < 6; ++q)
      value[q] = fmaf(weight[i], point[q], value[q]);
  }
}

__device__ __forceinline__ void
evaluate_field_table(const float *__restrict__ table, const Geometry &geometry,
                     double s, double theta, double phi, float value[4]) {
  int first_s, first_theta, first_phi;
  float weight_s[4], weight_theta[4], weight_phi[4];
  float derivative_theta[4], derivative_phi[4];
  cubic_table_location(s, 0, geometry, first_s, weight_s, nullptr);
  cubic_table_location(theta, 1, geometry, first_theta, weight_theta,
                       derivative_theta);
  cubic_table_location(phi, 2, geometry, first_phi, weight_phi, derivative_phi);
#pragma unroll
  for (int q = 0; q < 4; ++q)
    value[q] = 0.0f;

  const int ns = geometry.point_count[0];
  const int nt = geometry.point_count[1];
#pragma unroll
  for (int k = 0; k < 4; ++k) {
#pragma unroll
    for (int j = 0; j < 4; ++j) {
      const float weight_theta_phi = weight_theta[j] * weight_phi[k];
#pragma unroll
      for (int i = 0; i < 4; ++i) {
        const int point =
            (first_s + i) + ns * ((first_theta + j) + nt * (first_phi + k));
        const float bmod = table[2 * point];
        value[0] = fmaf(weight_s[i] * weight_theta_phi, bmod, value[0]);
        value[1] = fmaf(weight_s[i] * weight_theta_phi, table[2 * point + 1],
                        value[1]);
        value[2] = fmaf(weight_s[i] * derivative_theta[j] * weight_phi[k], bmod,
                        value[2]);
        value[3] = fmaf(weight_s[i] * weight_theta[j] * derivative_phi[k], bmod,
                        value[3]);
      }
    }
  }
}

__device__ __forceinline__ void
evaluate_boozer_table(const float *__restrict__ field_table,
                      const float *__restrict__ profile_table,
                      const Geometry &geometry, double s, double theta,
                      double phi, double &aphi, double &daphi, double &btheta,
                      double &dbtheta, double &bphi, double &dbphi,
                      double &bmod, double dbmod[3]) {
  float field_value[4], profile_value[6];
  evaluate_field_table(field_table, geometry, fabs(s), theta, phi, field_value);
  evaluate_profile_table(profile_table, geometry, fabs(s), profile_value);
  bmod = field_value[0];
  dbmod[0] = field_value[1];
  dbmod[1] = field_value[2];
  dbmod[2] = field_value[3];
  aphi = profile_value[0];
  daphi = profile_value[1];
  btheta = profile_value[2];
  dbtheta = profile_value[3];
  bphi = profile_value[4];
  dbphi = profile_value[5];
}

__device__ __noinline__ void
canonical_rhs(const float *__restrict__ field_table,
              const float *__restrict__ profile_table, const Geometry &geometry,
              const double y[4], double mu, double ro0, double dydt[4]) {
  const double radial = hypot(y[0], y[1]);
  const double theta = atan2(y[1], y[0]);
  double aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod[3];
  evaluate_boozer_table(field_table, profile_table, geometry, radial, theta,
                        y[2], aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod,
                        dbmod);

  const double inv_bmod = 1.0 / bmod;
  const double inv_bmod2 = inv_bmod * inv_bmod;
  const double htheta = btheta * inv_bmod;
  const double hphi = bphi * inv_bmod;
  const double dhtheta_r = dbtheta * inv_bmod - btheta * dbmod[0] * inv_bmod2;
  const double dhtheta_phi = -btheta * dbmod[2] * inv_bmod2;
  const double dhphi_r = dbphi * inv_bmod - bphi * dbmod[0] * inv_bmod2;
  const double dhphi_theta = -bphi * dbmod[1] * inv_bmod2;
  const double dhphi_phi = -bphi * dbmod[2] * inv_bmod2;
  const double inv_hphi = 1.0 / hphi;
  const double vparallel = (y[3] - aphi / ro0) * inv_hphi;
  const double dvparallel_r = -(daphi / ro0 + dhphi_r * vparallel) * inv_hphi;
  const double dvparallel_theta = -dhphi_theta * vparallel * inv_hphi;
  const double dvparallel_phi = -dhphi_phi * vparallel * inv_hphi;
  const double dh_r = vparallel * dvparallel_r + mu * dbmod[0];
  const double dh_theta = vparallel * dvparallel_theta + mu * dbmod[1];
  const double dh_phi = vparallel * dvparallel_phi + mu * dbmod[2];
  const double dptheta_r =
      dvparallel_r * htheta + vparallel * dhtheta_r + geometry.torflux / ro0;
  const double dptheta_phi = dvparallel_phi * htheta + vparallel * dhtheta_phi;
  const double inv_dptheta = 1.0 / dptheta_r;
  const double hprime = dh_r * inv_dptheta;

  const double radial_dot =
      -(dh_theta - htheta * inv_hphi * dh_phi) * inv_dptheta;
  const double theta_dot = hprime;
  dydt[0] =
      radial > 0.0 ? radial_dot * y[0] / radial - y[1] * theta_dot : radial_dot;
  dydt[1] = radial > 0.0 ? radial_dot * y[1] / radial + y[0] * theta_dot : 0.0;
  dydt[2] = (vparallel - hprime * htheta) * inv_hphi;
  dydt[3] = -(dh_phi - hprime * dptheta_phi);
}

__device__ __forceinline__ void
segment_limits(const float *__restrict__ field_table,
               const float *__restrict__ profile_table,
               const Geometry &geometry, const double y[4], double duration,
               double &hmax, double &momentum_atol_scale) {
  const double radial = hypot(y[0], y[1]);
  const double theta = atan2(y[1], y[0]);
  double aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod, dbmod[3];
  evaluate_boozer_table(field_table, profile_table, geometry, radial, theta,
                        y[2], aphi, daphi, btheta, dbtheta, bphi, dbphi, bmod,
                        dbmod);
  hmax = fmin(duration, 0.5 * kPi * fabs(bphi / bmod) / kSqrt2);
  momentum_atol_scale = kSqrt2 * fabs(bphi / bmod);
}

template <int Method> __device__ __forceinline__ int stage_count() {
  return Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE ? 7 : 6;
}

template <int Method>
__device__ __forceinline__ void stage_state(int stage, const double y[4],
                                            double h, const double k[7][4],
                                            double trial[4]) {
#pragma unroll
  for (int q = 0; q < 4; ++q)
    trial[q] = y[q];
  if constexpr (Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
    switch (stage) {
    case 1:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * k[0][q] / 5.0;
      break;
    case 2:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (3.0 * k[0][q] / 40.0 + 9.0 * k[1][q] / 40.0);
      break;
    case 3:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (44.0 * k[0][q] / 45.0 - 56.0 * k[1][q] / 15.0 +
                         32.0 * k[2][q] / 9.0);
      break;
    case 4:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] +=
            h * (19372.0 * k[0][q] / 6561.0 - 25360.0 * k[1][q] / 2187.0 +
                 64448.0 * k[2][q] / 6561.0 - 212.0 * k[3][q] / 729.0);
      break;
    case 5:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (9017.0 * k[0][q] / 3168.0 - 355.0 * k[1][q] / 33.0 +
                         46732.0 * k[2][q] / 5247.0 + 49.0 * k[3][q] / 176.0 -
                         5103.0 * k[4][q] / 18656.0);
      break;
    case 6:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (35.0 * k[0][q] / 384.0 + 500.0 * k[2][q] / 1113.0 +
                         125.0 * k[3][q] / 192.0 - 2187.0 * k[4][q] / 6784.0 +
                         11.0 * k[5][q] / 84.0);
      break;
    }
  } else {
    switch (stage) {
    case 1:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * k[0][q] / 5.0;
      break;
    case 2:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (3.0 * k[0][q] / 40.0 + 9.0 * k[1][q] / 40.0);
      break;
    case 3:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (3.0 * k[0][q] / 10.0 - 9.0 * k[1][q] / 10.0 +
                         6.0 * k[2][q] / 5.0);
      break;
    case 4:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] += h * (-11.0 * k[0][q] / 54.0 + 5.0 * k[1][q] / 2.0 -
                         70.0 * k[2][q] / 27.0 + 35.0 * k[3][q] / 27.0);
      break;
    case 5:
#pragma unroll
      for (int q = 0; q < 4; ++q)
        trial[q] +=
            h * (1631.0 * k[0][q] / 55296.0 + 175.0 * k[1][q] / 512.0 +
                 575.0 * k[2][q] / 13824.0 + 44275.0 * k[3][q] / 110592.0 +
                 253.0 * k[4][q] / 4096.0);
      break;
    }
  }
}

template <int Method>
__device__ __forceinline__ double
finish_attempt(const double y[4], double h, const double k[7][4],
               double trial[4], double tolerance, double momentum_atol_scale) {
  double error[4];
  if constexpr (Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
#pragma unroll
    for (int q = 0; q < 4; ++q) {
      trial[q] = y[q] + h * (35.0 * k[0][q] / 384.0 + 500.0 * k[2][q] / 1113.0 +
                             125.0 * k[3][q] / 192.0 -
                             2187.0 * k[4][q] / 6784.0 + 11.0 * k[5][q] / 84.0);
      error[q] = h * (71.0 * k[0][q] / 57600.0 - 71.0 * k[2][q] / 16695.0 +
                      71.0 * k[3][q] / 1920.0 - 17253.0 * k[4][q] / 339200.0 +
                      22.0 * k[5][q] / 525.0 - k[6][q] / 40.0);
    }
    double norm = 0.0;
#pragma unroll
    for (int q = 0; q < 4; ++q) {
      const double atol = q == 3 ? tolerance * momentum_atol_scale : tolerance;
      const double scale = atol + tolerance * (fabs(y[q]) + h * fabs(k[0][q]));
      norm = fmax(norm, fabs(error[q]) / scale);
    }
    return norm;
  } else {
#pragma unroll
    for (int q = 0; q < 4; ++q) {
      trial[q] =
          y[q] + h * (37.0 * k[0][q] / 378.0 + 250.0 * k[2][q] / 621.0 +
                      125.0 * k[3][q] / 594.0 + 512.0 * k[5][q] / 1771.0);
      const double fourth =
          y[q] + h * (2825.0 * k[0][q] / 27648.0 + 18575.0 * k[2][q] / 48384.0 +
                      13525.0 * k[3][q] / 55296.0 + 277.0 * k[4][q] / 14336.0 +
                      k[5][q] / 4.0);
      error[q] = trial[q] - fourth;
    }
    double sum = 0.0;
#pragma unroll
    for (int q = 0; q < 4; ++q) {
      const double atol = q == 3 ? tolerance * momentum_atol_scale : tolerance;
      const double scale = atol + tolerance * fmax(fabs(y[q]), fabs(trial[q]));
      const double ratio = error[q] / scale;
      sum = fma(ratio, ratio, sum);
    }
    return sqrt(sum / 4.0);
  }
}

template <int Method>
__device__ void trace_particle(const float *__restrict__ field_table,
                               const float *__restrict__ profile_table,
                               const Geometry &geometry,
                               const double initial_z[4], double mu, double ro0,
                               double final_z[4], double &lost_at,
                               int &particle_status, uint64_t &nfev) {
  double y[4] = {
      initial_z[0] * cos(initial_z[1]),
      initial_z[0] * sin(initial_z[1]),
      initial_z[2],
      initial_z[3],
  };
  double current_time = 0.0;
  lost_at = geometry.total_duration;
  particle_status = 0;
  nfev = 0;

  while (current_time < geometry.total_duration) {
    if (hypot(y[0], y[1]) >= 1.0) {
      lost_at = current_time;
      break;
    }
    const double duration =
        fmin(geometry.block_duration, geometry.total_duration - current_time);
    double hmax, momentum_atol_scale;
    segment_limits(field_table, profile_table, geometry, y, duration, hmax,
                   momentum_atol_scale);
    double h = fmin(duration, 1.0e-3 * hmax);
    double segment_time = 0.0;
    double previous_error = 1.0;
    bool first_step = true;
    bool after_reject = false;
    bool have_first_derivative = false;
    double first_derivative[4]{};

    for (int attempt = 0; attempt < kMaximumAttempts && segment_time < duration;
         ++attempt) {
      h = fmin(h, duration - segment_time);
      double k[7][4];
      if constexpr (Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
        if (!have_first_derivative) {
          canonical_rhs(field_table, profile_table, geometry, y, mu, ro0,
                        first_derivative);
          ++nfev;
          have_first_derivative = true;
        }
#pragma unroll
        for (int q = 0; q < 4; ++q)
          k[0][q] = first_derivative[q];
      } else {
        canonical_rhs(field_table, profile_table, geometry, y, mu, ro0, k[0]);
        ++nfev;
      }

      double trial[4];
      for (int stage = 1; stage < stage_count<Method>(); ++stage) {
        stage_state<Method>(stage, y, h, k, trial);
        canonical_rhs(field_table, profile_table, geometry, trial, mu, ro0,
                      k[stage]);
        ++nfev;
      }
      const double error = finish_attempt<Method>(
          y, h, k, trial, geometry.tolerance, momentum_atol_scale);
      if (!isfinite(error)) {
        particle_status = 1;
        break;
      }

      const double error_floor = fmax(error, 1.0e-300);
      double factor;
      if constexpr (Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
        factor = 0.9 * pow(error_floor, -1.0 / 3.0);
        factor = fmax(0.2, fmin(5.0, factor));
        if (error > 0.5 && error < 1.0)
          factor = 1.0;
      } else if (first_step || after_reject) {
        factor = 0.9 * pow(fmax(error, 1.0e-10), -1.0 / 5.0);
        factor = fmax(0.2, fmin(5.0, factor));
      } else {
        factor =
            0.9 * pow(fmax(error, 1.0e-10), -0.14) * pow(previous_error, 0.08);
        factor = fmax(0.2, fmin(5.0, factor));
      }
      const double hnew =
          fmax(geometry.minimum_timestep, fmin(hmax, h * factor));
      const bool at_minimum =
          geometry.minimum_timestep > 0.0 && h <= geometry.minimum_timestep;
      if (error <= 1.0 || at_minimum) {
        segment_time += h;
#pragma unroll
        for (int q = 0; q < 4; ++q)
          y[q] = trial[q];
        previous_error = fmax(error, 1.0e-10);
        first_step = false;
        after_reject = false;
        if constexpr (Method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
#pragma unroll
          for (int q = 0; q < 4; ++q)
            first_derivative[q] = k[6][q];
        }
        if (hypot(y[0], y[1]) >= 1.0) {
          lost_at = current_time + segment_time;
          break;
        }
      } else {
        after_reject = true;
      }
      h = hnew;
      if (attempt == kMaximumAttempts - 1)
        particle_status = 2;
    }
    if (particle_status != 0 || lost_at < geometry.total_duration)
      break;
    current_time += duration;
  }

  final_z[0] = hypot(y[0], y[1]);
  final_z[1] = atan2(y[1], y[0]);
  final_z[2] = y[2];
  final_z[3] = y[3];
}

template <int Method>
__global__ void
trace_kernel(const float *__restrict__ field_table,
             const float *__restrict__ profile_table, Geometry geometry,
             int particle_count, const double *__restrict__ initial_z,
             const double *__restrict__ mu, const double *__restrict__ ro0,
             double *__restrict__ final_z, double *__restrict__ loss_time,
             int *__restrict__ status, uint64_t *__restrict__ rhs_evaluations) {
  const int particle = blockIdx.x * blockDim.x + threadIdx.x;
  if (particle >= particle_count)
    return;
  trace_particle<Method>(field_table, profile_table, geometry,
                         initial_z + 4 * particle, mu[particle], ro0[particle],
                         final_z + 4 * particle, loss_time[particle],
                         status[particle], rhs_evaluations[particle]);
}

using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

int cuda_status(cudaError_t status) {
  return status == cudaSuccess ? 0 : 1000 + static_cast<int>(status);
}

} // namespace

extern "C" int simple_cuda_native_rk54(
    int method, int particle_count, const float *field_table,
    size_t field_table_count, const float *profile_table,
    size_t profile_table_count, const int point_count[3], const double x_min[3],
    const double inv_h_step[3], const double period[3],
    const double inv_period[3], double torflux, const double *initial_z,
    const double *mu, const double *ro0, double total_duration,
    double block_duration, double tolerance, double minimum_timestep,
    double *final_z, double *loss_time, int *status, uint64_t *rhs_evaluations,
    double profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_COUNT]) {
  if (profile_ms)
    std::fill(profile_ms, profile_ms + SIMPLE_CUDA_NATIVE_PROFILE_COUNT, 0.0);
  if (method != SIMPLE_CUDA_NATIVE_CASH_KARP &&
      method != SIMPLE_CUDA_NATIVE_DORMAND_PRINCE)
    return 1;
  if (particle_count <= 0 || !field_table || !profile_table || !point_count ||
      !x_min || !inv_h_step || !period || !inv_period || !initial_z || !mu ||
      !ro0 || !final_z || !loss_time || !status || !rhs_evaluations ||
      total_duration < 0.0 || block_duration <= 0.0 || tolerance <= 0.0 ||
      minimum_timestep < 0.0)
    return 2;

  size_t expected_field = 2;
  for (int d = 0; d < 3; ++d) {
    if (point_count[d] < 4 || inv_h_step[d] <= 0.0 ||
        (d > 0 && (period[d] <= 0.0 || inv_period[d] <= 0.0)))
      return 3;
    expected_field *= static_cast<size_t>(point_count[d]);
  }
  if (field_table_count != expected_field ||
      profile_table_count != 6 * static_cast<size_t>(point_count[0]))
    return 4;

  Geometry geometry{};
  for (int d = 0; d < 3; ++d) {
    geometry.point_count[d] = point_count[d];
    geometry.x_min[d] = x_min[d];
    geometry.inv_h_step[d] = inv_h_step[d];
    geometry.period[d] = period[d];
    geometry.inv_period[d] = inv_period[d];
  }
  geometry.torflux = torflux;
  geometry.total_duration = total_duration;
  geometry.block_duration = block_duration;
  geometry.tolerance = tolerance;
  geometry.minimum_timestep = minimum_timestep;

  const auto total_begin = Clock::now();
  const auto allocate_begin = total_begin;
  float *field_device = nullptr;
  float *profile_device = nullptr;
  double *initial_device = nullptr;
  double *mu_device = nullptr;
  double *ro0_device = nullptr;
  double *final_device = nullptr;
  double *loss_device = nullptr;
  int *status_device = nullptr;
  uint64_t *nfev_device = nullptr;
  cudaError_t error =
      cudaMalloc(&field_device, field_table_count * sizeof(float));
  if (error == cudaSuccess)
    error = cudaMalloc(&profile_device, profile_table_count * sizeof(float));
  if (error == cudaSuccess)
    error = cudaMalloc(&initial_device, 4 * particle_count * sizeof(double));
  if (error == cudaSuccess)
    error = cudaMalloc(&mu_device, particle_count * sizeof(double));
  if (error == cudaSuccess)
    error = cudaMalloc(&ro0_device, particle_count * sizeof(double));
  if (error == cudaSuccess)
    error = cudaMalloc(&final_device, 4 * particle_count * sizeof(double));
  if (error == cudaSuccess)
    error = cudaMalloc(&loss_device, particle_count * sizeof(double));
  if (error == cudaSuccess)
    error = cudaMalloc(&status_device, particle_count * sizeof(int));
  if (error == cudaSuccess)
    error = cudaMalloc(&nfev_device, particle_count * sizeof(uint64_t));
  const auto allocate_end = Clock::now();
  if (profile_ms)
    profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_ALLOCATE] =
        milliseconds(allocate_begin, allocate_end);
  if (error != cudaSuccess)
    goto cleanup;

  {
    const auto upload_begin = Clock::now();
    error =
        cudaMemcpy(field_device, field_table, field_table_count * sizeof(float),
                   cudaMemcpyHostToDevice);
    if (error == cudaSuccess)
      error = cudaMemcpy(profile_device, profile_table,
                         profile_table_count * sizeof(float),
                         cudaMemcpyHostToDevice);
    if (error == cudaSuccess)
      error = cudaMemcpy(initial_device, initial_z,
                         4 * particle_count * sizeof(double),
                         cudaMemcpyHostToDevice);
    if (error == cudaSuccess)
      error = cudaMemcpy(mu_device, mu, particle_count * sizeof(double),
                         cudaMemcpyHostToDevice);
    if (error == cudaSuccess)
      error = cudaMemcpy(ro0_device, ro0, particle_count * sizeof(double),
                         cudaMemcpyHostToDevice);
    const auto upload_end = Clock::now();
    if (profile_ms)
      profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_UPLOAD] =
          milliseconds(upload_begin, upload_end);
  }
  if (error != cudaSuccess)
    goto cleanup;

  {
    cudaEvent_t start{}, stop{};
    error = cudaEventCreate(&start);
    if (error == cudaSuccess)
      error = cudaEventCreate(&stop);
    if (error == cudaSuccess)
      error = cudaEventRecord(start);
    const int blocks = (particle_count + kThreads - 1) / kThreads;
    if (error == cudaSuccess) {
      if (method == SIMPLE_CUDA_NATIVE_DORMAND_PRINCE) {
        trace_kernel<SIMPLE_CUDA_NATIVE_DORMAND_PRINCE><<<blocks, kThreads>>>(
            field_device, profile_device, geometry, particle_count,
            initial_device, mu_device, ro0_device, final_device, loss_device,
            status_device, nfev_device);
      } else {
        trace_kernel<SIMPLE_CUDA_NATIVE_CASH_KARP><<<blocks, kThreads>>>(
            field_device, profile_device, geometry, particle_count,
            initial_device, mu_device, ro0_device, final_device, loss_device,
            status_device, nfev_device);
      }
      error = cudaGetLastError();
    }
    if (error == cudaSuccess)
      error = cudaEventRecord(stop);
    if (error == cudaSuccess)
      error = cudaEventSynchronize(stop);
    float kernel_ms = 0.0f;
    if (error == cudaSuccess)
      error = cudaEventElapsedTime(&kernel_ms, start, stop);
    if (profile_ms)
      profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_KERNEL] = kernel_ms;
    if (start)
      cudaEventDestroy(start);
    if (stop)
      cudaEventDestroy(stop);
  }
  if (error != cudaSuccess)
    goto cleanup;

  {
    const auto download_begin = Clock::now();
    error =
        cudaMemcpy(final_z, final_device, 4 * particle_count * sizeof(double),
                   cudaMemcpyDeviceToHost);
    if (error == cudaSuccess)
      error =
          cudaMemcpy(loss_time, loss_device, particle_count * sizeof(double),
                     cudaMemcpyDeviceToHost);
    if (error == cudaSuccess)
      error = cudaMemcpy(status, status_device, particle_count * sizeof(int),
                         cudaMemcpyDeviceToHost);
    if (error == cudaSuccess)
      error =
          cudaMemcpy(rhs_evaluations, nfev_device,
                     particle_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    const auto download_end = Clock::now();
    if (profile_ms)
      profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_DOWNLOAD] =
          milliseconds(download_begin, download_end);
  }

cleanup:
  cudaFree(nfev_device);
  cudaFree(status_device);
  cudaFree(loss_device);
  cudaFree(final_device);
  cudaFree(ro0_device);
  cudaFree(mu_device);
  cudaFree(initial_device);
  cudaFree(profile_device);
  cudaFree(field_device);
  if (profile_ms)
    profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_TOTAL] =
        milliseconds(total_begin, Clock::now());
  return cuda_status(error);
}

extern "C" const char *simple_cuda_native_error_string(int status) {
  switch (status) {
  case 0:
    return "success";
  case 1:
    return "unsupported integrator";
  case 2:
    return "invalid argument";
  case 3:
    return "invalid compact Boozer table geometry";
  case 4:
    return "compact Boozer table size mismatch";
  default:
    if (status >= 1000)
      return cudaGetErrorString(static_cast<cudaError_t>(status - 1000));
    return "unknown SIMPLE native CUDA error";
  }
}
