#include "simple_cuda_firm_rk54.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double angle_error(double actual, double expected) {
    return std::remainder(actual - expected, 2.0*pi);
}

int check_method(int method, bool cpu) {
    constexpr int particles = 16;
    constexpr int quantities = 13;
    constexpr int grid_points = 4;
    constexpr double mod_b = 2.0;
    constexpr double g = 3.0;
    constexpr double i = 0.4;
    constexpr double iota = 0.7;
    constexpr double total_speed = 1.0;
    constexpr double tmax = 0.02;
    const double ranges[9] = {
        0.0, 1.0, grid_points,
        0.0, pi, grid_points,
        0.0, 2.0*pi, grid_points
    };

    std::vector<double> field(64*quantities, 0.0);
    for (int point = 0; point < 64; ++point) {
        field[quantities*point] = mod_b;
        field[quantities*point + 4] = g;
        field[quantities*point + 6] = i;
        field[quantities*point + 8] = iota;
    }

    std::vector<double> initial_stz(3*particles);
    std::vector<double> initial_vparallel(particles);
    for (int particle = 0; particle < particles; ++particle) {
        initial_stz[3*particle] = 0.12 + 0.04*particle;
        initial_stz[3*particle + 1] = -2.5 + 0.31*particle;
        initial_stz[3*particle + 2] = -0.4 + 0.13*particle;
        initial_vparallel[particle] = 0.35 + 0.015*particle;
    }

    std::vector<double> output(7*particles);
    std::vector<unsigned long long> counters(3*particles);
    double profile[SIMPLE_CUDA_PROFILE_COUNT];
    int status;
    if (cpu) {
        status = simple_cpu_firm_rk54_landreman(
            method, particles, field.data(), field.size(), ranges,
            initial_stz.data(), initial_vparallel.data(), 4.0, 2.0,
            total_speed, 1.0, tmax, 1.0e-10, 1.0e-12, 0.02, 0.005,
            0.1, output.data(), counters.data(), profile);
    } else {
        status = simple_cuda_firm_rk54(
            method, particles, field.data(), field.size(), ranges,
            initial_stz.data(), initial_vparallel.data(), 4.0, 2.0,
            total_speed, 1.0, tmax, 1.0e-10, 1.0e-12, output.data(),
            counters.data(), profile);
    }
    if (status != 0) {
        std::fprintf(stderr, "%s RK method %d failed: %s\n",
                     cpu ? "CPU" : "CUDA", method,
                     simple_cuda_error_string(status));
        return 1;
    }

    double maximum_error = 0.0;
    for (int particle = 0; particle < particles; ++particle) {
        const double vparallel = initial_vparallel[particle];
        const double zeta_rate = vparallel*mod_b/(g + iota*i);
        const double theta_rate = iota*zeta_rate;
        const double expected_theta = initial_stz[3*particle + 1] +
                                      theta_rate*tmax;
        const double expected_zeta = initial_stz[3*particle + 2] +
                                     zeta_rate*tmax;
        maximum_error = std::max(maximum_error,
            std::fabs(output[7*particle] - tmax));
        maximum_error = std::max(maximum_error,
            std::fabs(output[7*particle + 1] - initial_stz[3*particle]));
        maximum_error = std::max(maximum_error,
            std::fabs(angle_error(output[7*particle + 2], expected_theta)));
        maximum_error = std::max(maximum_error,
            std::fabs(output[7*particle + 3] - expected_zeta));
        maximum_error = std::max(maximum_error,
            std::fabs(output[7*particle + 4] - vparallel));
        if (counters[3*particle] == 0 || counters[3*particle + 1] == 0) {
            std::fprintf(stderr, "%s RK method %d did not advance particle %d\n",
                         cpu ? "CPU" : "CUDA", method, particle);
            return 1;
        }
    }
    std::printf("backend=%s method=%d max analytic error=%.3e "
                "trace_ms=%.3f total_ms=%.3f\n",
                cpu ? "cpu" : "cuda", method, maximum_error,
                profile[SIMPLE_CUDA_PROFILE_KERNEL],
                profile[SIMPLE_CUDA_PROFILE_TOTAL]);
    if (maximum_error > 2.0e-9) {
        std::fprintf(stderr, "%s RK method %d exceeds analytic tolerance\n",
                     cpu ? "CPU" : "CUDA", method);
        return 1;
    }
    if (!(profile[SIMPLE_CUDA_PROFILE_TOTAL] >=
          profile[SIMPLE_CUDA_PROFILE_KERNEL])) {
        std::fprintf(stderr, "invalid CUDA profile timings\n");
        return 1;
    }
    return 0;
}

} // namespace

int main() {
    int failures = 0;
    const bool cpu_only = std::getenv("SIMPLE_CPU_ONLY") != nullptr;
    if (!cpu_only) {
        failures += check_method(SIMPLE_CUDA_DORMAND_PRINCE, false);
        failures += check_method(SIMPLE_CUDA_CASH_KARP, false);
    }
    failures += check_method(SIMPLE_CUDA_DORMAND_PRINCE, true);
    failures += check_method(SIMPLE_CUDA_CASH_KARP, true);
    if (cpu_only) return failures == 0 ? 0 : 1;

    const double ranges[9] = {0.0, 1.0, 4.0, 0.0, pi, 4.0,
                              0.0, 2.0*pi, 4.0};
    double field[64*13]{};
    double stz[3]{};
    double vparallel[1]{};
    double output[7]{};
    unsigned long long counters[3]{};
    double profile[SIMPLE_CUDA_PROFILE_COUNT]{};
    const int invalid = simple_cuda_firm_rk54(
        SIMPLE_CUDA_DORMAND_PRINCE, 1, field, 64*13 - 1, ranges, stz,
        vparallel, 1.0, 1.0, 1.0, 1.0, 0.1, 1.0e-6, 0.0, output,
        counters, profile);
    if (invalid != 4) {
        std::fprintf(stderr, "metagrid-size validation returned %d, expected 4\n",
                     invalid);
        ++failures;
    }

    for (int point = 0; point < 64; ++point) {
        field[13*point] = 2.0;
        field[13*point + 4] = 3.0;
        field[13*point + 6] = 0.4;
        field[13*point + 8] = 0.7;
    }
    stz[0] = 0.5;
    stz[1] = 0.2;
    vparallel[0] = 0.4;
    const int segmented = simple_cuda_firm_rk54_landreman(
        SIMPLE_CUDA_DORMAND_PRINCE, 1, field, 64*13, ranges, stz,
        vparallel, 4.0, 2.0, 1.0, 1.0, 0.02, 1.0e-10, 1.0e-12,
        0.02, 0.005, 0.1, output, counters, profile);
    if (segmented != 0 || std::fabs(output[0] - 0.02) > 1.0e-15 ||
        counters[1] == 0) {
        std::fprintf(stderr,
            "Landreman segmented path failed: status=%d time=%.17g accepted=%llu\n",
            segmented, output[0], counters[1]);
        ++failures;
    }
    return failures == 0 ? 0 : 1;
}
