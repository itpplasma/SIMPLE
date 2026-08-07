#ifndef SIMPLE_CUDA_FIRM_RK54_H
#define SIMPLE_CUDA_FIRM_RK54_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
    SIMPLE_CUDA_CASH_KARP = 1,
    SIMPLE_CUDA_DORMAND_PRINCE = 2
};

enum {
    SIMPLE_CUDA_PROFILE_ALLOCATE = 0,
    SIMPLE_CUDA_PROFILE_UPLOAD = 1,
    SIMPLE_CUDA_PROFILE_KERNEL = 2,
    SIMPLE_CUDA_PROFILE_DOWNLOAD = 3,
    SIMPLE_CUDA_PROFILE_METRIC = 4,
    SIMPLE_CUDA_PROFILE_TOTAL = 5,
    SIMPLE_CUDA_PROFILE_COUNT = 6
};

/*
 * Trace guiding centers in the same pseudocartesian Boozer state and the same
 * 13-quantity metagrid layout used by FIRM3D's Landreman GPU calculation.
 *
 * ranges contains three consecutive [minimum, maximum, number_of_points]
 * triples. quad_points is laid out as
 *   [cell_s][cell_theta][cell_zeta][4][4][4][13],
 * with the final index contiguous. initial_stz and final_stzv are particle-
 * major. final_stzv holds [loss_time, s, theta, zeta, v_parallel, dt_min,
 * dt_max]. counters holds [rhs_evaluations, accepted_steps, rejected_steps]
 * per particle. profile_ms contains the six timings indexed above; TOTAL also
 * includes cleanup and other host overhead, so it is not an exclusive phase.
 */
int simple_cuda_firm_rk54(
    int method,
    int particle_count,
    const double *quad_points,
    size_t quad_point_count,
    const double ranges[9],
    const double *initial_stz,
    const double *initial_vparallel,
    double mass,
    double charge,
    double total_speed,
    double psi0,
    double tmax,
    double tolerance,
    double minimum_timestep,
    double *final_stzv,
    unsigned long long *counters,
    double profile_ms[SIMPLE_CUDA_PROFILE_COUNT]);

/* Landreman objective variant: keep the field resident, trace in t_block
 * segments, and stop after the exp(-t/tau)-weighted loss exceeds maxloss.
 * This reproduces FIRM3D's host-side early-stopping protocol. */
int simple_cuda_firm_rk54_landreman(
    int method,
    int particle_count,
    const double *quad_points,
    size_t quad_point_count,
    const double ranges[9],
    const double *initial_stz,
    const double *initial_vparallel,
    double mass,
    double charge,
    double total_speed,
    double psi0,
    double tmax,
    double tolerance,
    double minimum_timestep,
    double maxloss,
    double t_block,
    double tau,
    double *final_stzv,
    unsigned long long *counters,
    double profile_ms[SIMPLE_CUDA_PROFILE_COUNT]);

const char *simple_cuda_error_string(int status);

#ifdef __cplusplus
}
#endif

#endif
