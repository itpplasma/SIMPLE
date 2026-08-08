#include "simple_cuda_native_rk54.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {

bool close(double actual, double expected, double tolerance = 1.0e-7) {
  return std::abs(actual - expected) <= tolerance;
}

int exercise_method(int method) {
  constexpr int particle_count = 3;
  const int point_count[3] = {4, 4, 4};
  const double x_min[3] = {0.0, 0.0, 0.0};
  const double inv_h_step[3] = {1.0, 1.0, 1.0};
  const double period[3] = {0.0, 3.0, 3.0};
  const double inv_period[3] = {0.0, 1.0 / 3.0, 1.0 / 3.0};
  std::vector<float> field_table(2 * 4 * 4 * 4, 0.0f);
  for (std::size_t point = 0; point < field_table.size() / 2; ++point)
    field_table[2 * point] = 1.0f;
  std::vector<float> profile_table(6 * 4, 0.0f);
  for (int radial = 0; radial < 4; ++radial)
    profile_table[6 * radial + 4] = 1.0f;

  const std::array<double, 4 * particle_count> initial_z = {
      0.2, 0.3, 0.4, 0.5, 0.4, -0.2, 1.1, -0.7, 0.6, 0.8, -0.5, 1.25,
  };
  const std::array<double, particle_count> mu = {0.1, 0.2, 0.3};
  const std::array<double, particle_count> ro0 = {2.0, 3.0, 4.0};
  std::array<double, 4 * particle_count> final_z{};
  std::array<double, particle_count> loss_time{};
  std::array<int, particle_count> status{};
  std::array<uint64_t, particle_count> nfev{};
  uint64_t warp_rhs_slots{};
  double observed_duration{};
  double profile[SIMPLE_CUDA_NATIVE_PROFILE_COUNT]{};
  constexpr double duration = 0.02;

  const int result = simple_cuda_native_rk54(
      method, particle_count, field_table.data(), field_table.size(),
      profile_table.data(), profile_table.size(), point_count, x_min,
      inv_h_step, period, inv_period, 1.0, initial_z.data(), mu.data(),
      ro0.data(), duration, 0.005, 1.0e-10, 1.0e-14, 1.0, 1.0,
      &observed_duration, final_z.data(), loss_time.data(), status.data(),
      nfev.data(), &warp_rhs_slots, profile);
  if (result != 0) {
    std::fprintf(stderr, "%s\n", simple_cuda_native_error_string(result));
    return 1;
  }
  if (!close(observed_duration, duration, 1.0e-14)) {
    std::fprintf(stderr, "unexpected observed duration %.17g\n",
                 observed_duration);
    return 2;
  }
  for (int particle = 0; particle < particle_count; ++particle) {
    const int offset = 4 * particle;
    if (!close(final_z[offset], initial_z[offset]) ||
        !close(final_z[offset + 1], initial_z[offset + 1]) ||
        !close(final_z[offset + 2],
               initial_z[offset + 2] + initial_z[offset + 3] * duration) ||
        !close(final_z[offset + 3], initial_z[offset + 3]) ||
        !close(loss_time[particle], duration) || status[particle] != 0 ||
        nfev[particle] == 0) {
      std::fprintf(stderr,
                   "closed-form mismatch method %d particle %d: "
                   "actual [%.17g, %.17g, %.17g, %.17g], expected phi %.17g, "
                   "loss %.17g, status %d, nfev %llu\n",
                   method, particle, final_z[offset], final_z[offset + 1],
                   final_z[offset + 2], final_z[offset + 3],
                   initial_z[offset + 2] + initial_z[offset + 3] * duration,
                   loss_time[particle], status[particle],
                   static_cast<unsigned long long>(nfev[particle]));
      return 2;
    }
  }
  if (!(profile[SIMPLE_CUDA_NATIVE_PROFILE_KERNEL] > 0.0)) {
    std::fprintf(stderr, "missing CUDA kernel timing\n");
    return 3;
  }

  const int invalid = simple_cuda_native_rk54(
      method, particle_count, field_table.data(), field_table.size() - 1,
      profile_table.data(), profile_table.size(), point_count, x_min,
      inv_h_step, period, inv_period, 1.0, initial_z.data(), mu.data(),
      ro0.data(), duration, 0.005, 1.0e-10, 1.0e-14, 1.0, 1.0,
      &observed_duration, final_z.data(), loss_time.data(), status.data(),
      nfev.data(), &warp_rhs_slots, profile);
  if (invalid != 4) {
    std::fprintf(stderr,
                 "table-size validation did not reject malformed input\n");
    return 4;
  }
  return 0;
}

int exercise_early_stop(int method) {
  constexpr int particle_count = 2;
  const int point_count[3] = {4, 4, 4};
  const double x_min[3] = {0.0, 0.0, 0.0};
  const double inv_h_step[3] = {1.0, 1.0, 1.0};
  const double period[3] = {0.0, 3.0, 3.0};
  const double inv_period[3] = {0.0, 1.0 / 3.0, 1.0 / 3.0};
  std::vector<float> field_table(2 * 4 * 4 * 4, 0.0f);
  for (std::size_t point = 0; point < field_table.size() / 2; ++point)
    field_table[2 * point] = 1.0f;
  std::vector<float> profile_table(6 * 4, 0.0f);
  for (int radial = 0; radial < 4; ++radial)
    profile_table[6 * radial + 4] = 1.0f;

  const std::array<double, 4 * particle_count> initial_z = {
      1.1, 0.0, 0.0, 0.5, 0.2, 0.3, 0.4, -0.7,
  };
  const std::array<double, particle_count> mu = {0.1, 0.2};
  const std::array<double, particle_count> ro0 = {2.0, 3.0};
  std::array<double, 4 * particle_count> final_z{};
  std::array<double, particle_count> loss_time{};
  std::array<int, particle_count> status{};
  std::array<uint64_t, particle_count> nfev{};
  uint64_t warp_rhs_slots{};
  double observed_duration{};
  double profile[SIMPLE_CUDA_NATIVE_PROFILE_COUNT]{};
  constexpr double duration = 0.02;
  constexpr double block_duration = 0.005;

  const int result = simple_cuda_native_rk54(
      method, particle_count, field_table.data(), field_table.size(),
      profile_table.data(), profile_table.size(), point_count, x_min,
      inv_h_step, period, inv_period, 1.0, initial_z.data(), mu.data(),
      ro0.data(), duration, block_duration, 1.0e-10, 1.0e-14, 1.0, 0.0,
      &observed_duration, final_z.data(), loss_time.data(), status.data(),
      nfev.data(), &warp_rhs_slots, profile);
  const double expected_survivor_phi =
      initial_z[6] + initial_z[7] * block_duration;
  if (result != 0 || !close(observed_duration, block_duration, 1.0e-14) ||
      !close(loss_time[0], 0.0, 1.0e-14) ||
      !close(loss_time[1], duration, 1.0e-14) || status[0] != 0 ||
      status[1] != 0 || !close(final_z[6], expected_survivor_phi)) {
    std::fprintf(stderr,
                 "early-stop mismatch method %d: result %d, observed %.17g, "
                 "losses [%.17g, %.17g], survivor phi %.17g\n",
                 method, result, observed_duration, loss_time[0], loss_time[1],
                 final_z[6]);
    return 5;
  }
  return 0;
}

} // namespace

int main() {
  for (const int method :
       {SIMPLE_CUDA_NATIVE_DORMAND_PRINCE, SIMPLE_CUDA_NATIVE_CASH_KARP}) {
    const int result = exercise_method(method);
    if (result != 0)
      return result;
    const int early_stop_result = exercise_early_stop(method);
    if (early_stop_result != 0)
      return early_stop_result;
  }
  std::printf("native CUDA RK54 closed-form oracle passed\n");
  return 0;
}
