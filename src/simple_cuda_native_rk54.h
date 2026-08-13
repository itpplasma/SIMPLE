#ifndef SIMPLE_CUDA_NATIVE_RK54_H
#define SIMPLE_CUDA_NATIVE_RK54_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  SIMPLE_CUDA_NATIVE_CASH_KARP = 1,
  /* Public DOPRI production dispatch uses the tuned PI controller. */
  SIMPLE_CUDA_NATIVE_DORMAND_PRINCE_TUNED = 2,
  /* Keep the old spelling and numeric ABI for existing callers. */
  SIMPLE_CUDA_NATIVE_DORMAND_PRINCE = SIMPLE_CUDA_NATIVE_DORMAND_PRINCE_TUNED,
};

enum {
  SIMPLE_CUDA_NATIVE_PROFILE_ALLOCATE = 0,
  SIMPLE_CUDA_NATIVE_PROFILE_UPLOAD = 1,
  SIMPLE_CUDA_NATIVE_PROFILE_KERNEL = 2,
  SIMPLE_CUDA_NATIVE_PROFILE_DOWNLOAD = 3,
  SIMPLE_CUDA_NATIVE_PROFILE_TOTAL = 4,
  SIMPLE_CUDA_NATIVE_PROFILE_COUNT = 5,
};

int simple_cuda_native_rk54(
    int method, int particle_count, const float *field_table,
    size_t field_table_count, const float *profile_table,
    size_t profile_table_count, const int point_count[3], const double x_min[3],
    const double inv_h_step[3], const double period[3],
    const double inv_period[3], double torflux, const double *initial_z,
    const double *mu, const double *ro0, double total_duration,
    double block_duration, double tolerance, double minimum_timestep,
    double loss_decay_rate, double maxloss, double *observed_duration,
    double *final_z, double *loss_time, int *status, uint64_t *rhs_evaluations,
    uint64_t *warp_rhs_slots,
    double profile_ms[SIMPLE_CUDA_NATIVE_PROFILE_COUNT]);

const char *simple_cuda_native_error_string(int status);

#ifdef __cplusplus
}
#endif

#endif
