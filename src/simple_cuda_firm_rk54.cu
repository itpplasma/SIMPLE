#include "simple_cuda_firm_rk54.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace {

constexpr int kStoredQuantities = 13;
#ifndef SIMPLE_CUDA_THREADS
#define SIMPLE_CUDA_THREADS 32
#endif
constexpr int kThreads = SIMPLE_CUDA_THREADS;
static_assert(kThreads > 0 && kThreads <= 1024,
              "SIMPLE_CUDA_THREADS must be between 1 and 1024");
constexpr int kMaximumAttempts = 100000000;
constexpr int kMaximumStages = 13;
constexpr double kPi = 3.141592653589793238462643383279502884;

template<int Method>
__host__ __device__ constexpr bool is_dopri5() {
    return Method == SIMPLE_CUDA_DORMAND_PRINCE ||
           Method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED;
}

template<int Method>
__host__ __device__ constexpr bool is_dop853() {
    return Method == SIMPLE_CUDA_DOP853;
}

template<int Method>
__host__ __device__ constexpr bool uses_dopri_field_layout() {
    return is_dopri5<Method>() || is_dop853<Method>();
}

bool supported_method(int method) {
    return method == SIMPLE_CUDA_CASH_KARP ||
           method == SIMPLE_CUDA_DORMAND_PRINCE ||
           method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED ||
           method == SIMPLE_CUDA_DOP853;
}

struct Parameters {
    double minimum[3];
    double period_or_maximum[3];
    double step[3];
    int point_count[3];
    int cell_count[3];
    double mass;
    double charge;
    double total_speed;
    double psi0;
    double tmax;
    double tolerance;
    double minimum_timestep;
};

__host__ __device__ __forceinline__ double wrap_periodic(double x,
                                                          double minimum,
                                                          double period) {
    double wrapped = fmod(x - minimum, period);
    if (wrapped < 0.0) wrapped += period;
    return minimum + wrapped;
}

__host__ __device__ __forceinline__ void lagrange_weights(double x,
                                                           double w[4]) {
    w[0] = (1.0 - x)*(2.0 - x)*(3.0 - x)/6.0;
    w[1] = x*(2.0 - x)*(3.0 - x)/2.0;
    w[2] = x*(x - 1.0)*(3.0 - x)/2.0;
    w[3] = x*(x - 1.0)*(x - 2.0)/6.0;
}

template<int Method>
__host__ __device__ __noinline__ void interpolate13(
    const double *__restrict__ data, const Parameters &p,
    const double state[4], const double s, double theta,
    double *output, bool &reflected) {
    constexpr int quantity_count = uses_dopri_field_layout<Method>() ? 12 : 13;
    // The stored theta metagrid is the stellarator-symmetric half-domain
    // [0, pi], but the physical angle must first be wrapped over 2*pi and only
    // then reflected. Wrapping directly by the metagrid extent changes the
    // sign of odd quantities and is not equivalent.
    if (theta < 0.0) theta += 2.0*kPi;
    double zeta = wrap_periodic(state[2], p.minimum[2],
                                p.period_or_maximum[2] - p.minimum[2]);

    reflected = theta > kPi;
    if (reflected) {
        theta = 2.0*kPi - theta;
        zeta = p.period_or_maximum[2] - zeta;
        if (zeta >= p.period_or_maximum[2]) zeta = p.minimum[2];
    }

    const double coordinate[3] = {s, theta, zeta};
    int grid_index[3];
    double relative[3];
#pragma unroll
    for (int d = 0; d < 3; ++d) {
        int first = 3*(static_cast<int>((coordinate[d] - p.minimum[d])/
                                        p.step[d])/3);
        first = first < 0 ? 0 : first;
        first = first > p.point_count[d] - 4 ? p.point_count[d] - 4 : first;
        grid_index[d] = first/3;
        relative[d] = (coordinate[d] - p.minimum[d] - first*p.step[d])/
                      p.step[d];
    }

    double w0[4], w1[4], w2[4];
    lagrange_weights(relative[0], w0);
    lagrange_weights(relative[1], w1);
    lagrange_weights(relative[2], w2);

    const std::size_t cell =
        (static_cast<std::size_t>(grid_index[0])*p.cell_count[1] +
         grid_index[1])*p.cell_count[2] + grid_index[2];
    const double *cell_data = data + kStoredQuantities*64*cell;

#pragma unroll
    for (int q = 0; q < quantity_count; ++q) output[q] = 0.0;
#pragma unroll
    for (int i = 0; i < 4; ++i) {
#pragma unroll
        for (int j = 0; j < 4; ++j) {
            const double wij = w0[i]*w1[j];
#pragma unroll
            for (int k = 0; k < 4; ++k) {
                const double weight = wij*w2[k];
                const double *point =
                    cell_data + kStoredQuantities*(16*i + 4*j + k);
                if constexpr (uses_dopri_field_layout<Method>()) {
#pragma unroll
                    for (int q = 0; q < 9; ++q)
                        output[q] = fma(point[q], weight, output[q]);
                    // d(iota)/d(psi), stored at point[9], does not enter the
                    // general guiding-center equations.
                    output[9] = fma(point[10], weight, output[9]);
                    output[10] = fma(point[11], weight, output[10]);
                    output[11] = fma(point[12], weight, output[11]);
                } else {
#pragma unroll
                    for (int q = 0; q < kStoredQuantities; ++q)
                        output[q] = fma(point[q], weight, output[q]);
                }
            }
        }
    }
}

template<int Method>
__host__ __device__ __forceinline__ void rhs_from_interpolated(
    const Parameters &p, const double mu, const double state[4],
    const double value[], bool reflected, double derivative[4]) {
    const double x = state[0];
    const double y = state[1];
    const double s = hypot(x, y);
    const double vparallel = state[3];

    const double mod_b = value[0];
    const double db_dpsi = value[1]/p.psi0;
    double db_dtheta = value[2];
    double db_dzeta = value[3];
    const double g = value[4];
    const double dg_dpsi = value[5]/p.psi0;
    const double i = value[6];
    const double di_dpsi = value[7]/p.psi0;
    const double iota = value[8];
    double k;
    double dk_dtheta;
    double dk_dzeta;
    if constexpr (uses_dopri_field_layout<Method>()) {
        k = value[9];
        dk_dtheta = value[10];
        dk_dzeta = value[11];
    } else {
        k = value[10];
        dk_dtheta = value[11];
        dk_dzeta = value[12];
    }
    if (reflected) {
        db_dtheta = -db_dtheta;
        db_dzeta = -db_dzeta;
        k = -k;
    }

    const double c = -p.mass*vparallel*dk_dzeta/mod_b - p.charge*iota +
                     p.mass*vparallel*dg_dpsi/mod_b;
    const double f = -p.mass*vparallel*dk_dtheta/mod_b + p.charge +
                     p.mass*vparallel*di_dpsi/mod_b;
    const double d = (f*g - c*i)/iota;
    const double energy_factor = p.mass*(vparallel*vparallel/mod_b + mu);
    const double sdot = (i*db_dzeta - g*db_dtheta)*energy_factor/
                        (iota*d*p.psi0);
    const double thetadot = ((g*db_dpsi - k*db_dzeta)*energy_factor -
                             c*vparallel*mod_b)/(iota*d);
    const double zetadot = (f*vparallel*mod_b -
                            (db_dpsi*i - db_dtheta*k)*energy_factor)/(iota*d);
    const double vpardot = (c*db_dtheta - f*db_dzeta)*mu*mod_b/(iota*d);

    const double radial_x = s > 0.0 ? x/s : 1.0;
    const double radial_y = s > 0.0 ? y/s : 0.0;
    derivative[0] = sdot*radial_x - y*thetadot;
    derivative[1] = sdot*radial_y + x*thetadot;
    derivative[2] = zetadot;
    derivative[3] = vpardot;
}

template<int Method>
__host__ __device__ __noinline__ void right_hand_side(
    const double *__restrict__ field, const Parameters &p, const double mu,
    const double state[4], double derivative[4], double *mod_b_out = nullptr,
    double *g_out = nullptr) {
    const double x = state[0];
    const double y = state[1];
    const double s = hypot(x, y);
    const double theta = atan2(y, x);
    constexpr int quantity_count = uses_dopri_field_layout<Method>() ? 12 : 13;
    double value[quantity_count];
    bool reflected;
    interpolate13<Method>(field, p, state, s, theta, value, reflected);
    rhs_from_interpolated<Method>(p, mu, state, value, reflected, derivative);
    const double mod_b = value[0];
    const double g = value[4];
    if (mod_b_out) *mod_b_out = mod_b;
    if (g_out) *g_out = g;
}

template<int Method>
__host__ __device__ __noinline__ void initialize_particle(
    const double *__restrict__ field, const Parameters &p, const double state[4],
    double &mu, double &hmax, double derivative[4]) {
    const double x = state[0];
    const double y = state[1];
    const double s = hypot(x, y);
    const double theta = atan2(y, x);
    constexpr int quantity_count = uses_dopri_field_layout<Method>() ? 12 : 13;
    double value[quantity_count];
    bool reflected;
    interpolate13<Method>(field, p, state, s, theta, value, reflected);
    const double mod_b = value[0];
    mu = (p.total_speed*p.total_speed - state[3]*state[3])/(2.0*mod_b);
    hmax = fmin(p.tmax, (value[4]/mod_b)*(0.5*kPi)/p.total_speed);
    rhs_from_interpolated<Method>(p, mu, state, value, reflected, derivative);
}

template<int Method>
__host__ __device__ __forceinline__ int stage_count() {
    if constexpr (is_dop853<Method>()) return 12;
    return is_dopri5<Method>() ? 7 : 6;
}

template<int Method>
__host__ __device__ __forceinline__ void stage_state(
    int stage, const double y[4], double h,
    double k[kMaximumStages][4], double trial[4]) {
#pragma unroll
    for (int q = 0; q < 4; ++q) trial[q] = y[q];
    if constexpr (is_dop853<Method>()) {
        switch (stage) {
        case 1:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*0.05260015195876773*k[0][q];
            break;
        case 2:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.01972505698453790*k[0][q] +
                               0.05917517095361370*k[1][q]);
            break;
        case 3:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.02958758547680685*k[0][q] +
                               0.08876275643042055*k[2][q]);
            break;
        case 4:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.24136513415926669*k[0][q] -
                               0.88454947932828609*k[2][q] +
                               0.92483400326179199*k[3][q]);
            break;
        case 5:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.03703703703703704*k[0][q] +
                               0.17082860872947386*k[3][q] +
                               0.12546768756682242*k[4][q]);
            break;
        case 6:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.03710937500000000*k[0][q] +
                               0.17025221101954405*k[3][q] +
                               0.06021653898045596*k[4][q] -
                               0.01757812500000000*k[5][q]);
            break;
        case 7:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.03709200011850479*k[0][q] +
                               0.17038392571223998*k[3][q] +
                               0.10726203044637328*k[4][q] -
                               0.01531943774862441*k[5][q] +
                               0.00827378916381402*k[6][q]);
            break;
        case 8:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.62411095871607569*k[0][q] -
                               3.36089262944694142*k[3][q] -
                               0.86821934684172597*k[4][q] +
                               27.59209969944671004*k[5][q] +
                               20.15406755047789014*k[6][q] -
                               43.48988418106996111*k[7][q]);
            break;
        case 9:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(0.47766253643826434*k[0][q] -
                               2.48811461997166770*k[3][q] -
                               0.59029082683684297*k[4][q] +
                               21.23005144818119311*k[5][q] +
                               15.27923363288242315*k[6][q] -
                               33.28821096898486294*k[7][q] -
                               0.02033120170850863*k[8][q]);
            break;
        case 10:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(-0.93714243008598730*k[0][q] +
                               5.18637242884406380*k[3][q] +
                               1.09143734899672950*k[4][q] -
                               8.14978701074692680*k[5][q] -
                               18.520065659996959*k[6][q] +
                               22.739487099350505*k[7][q] +
                               2.49360555267965230*k[8][q] -
                               3.04676447189821960*k[9][q]);
            break;
        case 11:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(2.27331014751653800*k[0][q] -
                               10.534495466737249*k[3][q] -
                               2.00087205822486250*k[4][q] -
                               17.958931863118799*k[5][q] +
                               27.948884529419960*k[6][q] -
                               2.85899827713502350*k[7][q] -
                               8.87285693353062930*k[8][q] +
                               12.360567175794303*k[9][q] +
                               0.64339274601576357*k[10][q]);
            break;
        }
    } else if constexpr (is_dopri5<Method>()) {
        switch (stage) {
        case 1:
#pragma unroll
            for (int q = 0; q < 4; ++q) trial[q] += h*(k[0][q]/5.0);
            break;
        case 2:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(3.0*k[0][q]/40.0 + 9.0*k[1][q]/40.0);
            break;
        case 3:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(44.0*k[0][q]/45.0 - 56.0*k[1][q]/15.0 +
                               32.0*k[2][q]/9.0);
            break;
        case 4:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(19372.0*k[0][q]/6561.0 -
                               25360.0*k[1][q]/2187.0 +
                               64448.0*k[2][q]/6561.0 -
                               212.0*k[3][q]/729.0);
            break;
        case 5:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(9017.0*k[0][q]/3168.0 -
                               355.0*k[1][q]/33.0 +
                               46732.0*k[2][q]/5247.0 +
                               49.0*k[3][q]/176.0 -
                               5103.0*k[4][q]/18656.0);
            break;
        case 6:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(35.0*k[0][q]/384.0 +
                               500.0*k[2][q]/1113.0 +
                               125.0*k[3][q]/192.0 -
                               2187.0*k[4][q]/6784.0 +
                               11.0*k[5][q]/84.0);
            break;
        }
    } else {
        switch (stage) {
        case 1:
#pragma unroll
            for (int q = 0; q < 4; ++q) trial[q] += h*k[0][q]/5.0;
            break;
        case 2:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(3.0*k[0][q]/40.0 + 9.0*k[1][q]/40.0);
            break;
        case 3:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(3.0*k[0][q]/10.0 - 9.0*k[1][q]/10.0 +
                               6.0*k[2][q]/5.0);
            break;
        case 4:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(-11.0*k[0][q]/54.0 + 5.0*k[1][q]/2.0 -
                               70.0*k[2][q]/27.0 + 35.0*k[3][q]/27.0);
            break;
        case 5:
#pragma unroll
            for (int q = 0; q < 4; ++q)
                trial[q] += h*(1631.0*k[0][q]/55296.0 +
                               175.0*k[1][q]/512.0 +
                               575.0*k[2][q]/13824.0 +
                               44275.0*k[3][q]/110592.0 +
                               253.0*k[4][q]/4096.0);
            break;
        }
    }
}

template<int Method>
__host__ __device__ __forceinline__ double finish_attempt(
    const double y[4], double h, double k[kMaximumStages][4], double trial[4],
    double tolerance, double total_speed) {
    double error[4];
    if constexpr (is_dop853<Method>()) {
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            trial[q] = y[q] + h*(0.05429373411656876*k[0][q] +
                4.45031289275240900*k[5][q] +
                1.89151789931450030*k[6][q] -
                5.80120396001058500*k[7][q] +
                0.31116436695781990*k[8][q] -
                0.15216094966251610*k[9][q] +
                0.20136540080403034*k[10][q] +
                0.04471061572777259*k[11][q]);
            const double error5 = h*(0.01312004499419488*k[0][q] -
                1.22515644637620440*k[5][q] -
                0.49575894965725020*k[6][q] +
                1.66437718245498640*k[7][q] -
                0.35032884874997366*k[8][q] +
                0.33417911871301750*k[9][q] +
                0.08192320648511571*k[10][q] -
                0.02235530786388629*k[11][q]);
            const double atol = q == 3 ? tolerance*total_speed : tolerance;
            const double scale = atol + tolerance*fmax(fabs(y[q]),
                                                       fabs(trial[q]));
            error[q] = error5/scale;
        }
        double error5_norm_sq = 0.0;
        double error3_norm_sq = 0.0;
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            const double error5 = error[q];
            const double error3 = h*( -0.18980075407240762*k[0][q] +
                4.45031289275240900*k[5][q] +
                1.89151789931450030*k[6][q] -
                5.80120396001058500*k[7][q] -
                0.42268232132379190*k[8][q] -
                0.15216094966251610*k[9][q] +
                0.20136540080403034*k[10][q] +
                0.02265179219836082*k[11][q]);
            const double component_atol = q == 3 ? tolerance*total_speed :
                                           tolerance;
            const double scale = component_atol +
                tolerance*fmax(fabs(y[q]), fabs(trial[q]));
            const double scaled3 = error3/scale;
            error5_norm_sq = fma(error5, error5, error5_norm_sq);
            error3_norm_sq = fma(scaled3, scaled3, error3_norm_sq);
        }
        if (error5_norm_sq == 0.0 && error3_norm_sq == 0.0) return 0.0;
        return error5_norm_sq/
               sqrt((error5_norm_sq + 0.01*error3_norm_sq)*4.0);
    } else if constexpr (is_dopri5<Method>()) {
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            trial[q] = y[q] + h*(35.0*k[0][q]/384.0 +
                500.0*k[2][q]/1113.0 + 125.0*k[3][q]/192.0 -
                2187.0*k[4][q]/6784.0 + 11.0*k[5][q]/84.0);
            error[q] = h*(71.0*k[0][q]/57600.0 - 71.0*k[2][q]/16695.0 +
                71.0*k[3][q]/1920.0 - 17253.0*k[4][q]/339200.0 +
                22.0*k[5][q]/525.0 - k[6][q]/40.0);
        }
        double norm = 0.0;
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            const double atol = q == 3 ? tolerance*total_speed : tolerance;
            const double scale = atol + tolerance*(fabs(y[q]) + h*fabs(k[0][q]));
            norm = fmax(norm, fabs(error[q])/scale);
        }
        return norm;
    } else {
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            trial[q] = y[q] + h*(37.0*k[0][q]/378.0 +
                250.0*k[2][q]/621.0 + 125.0*k[3][q]/594.0 +
                512.0*k[5][q]/1771.0);
            const double fourth = y[q] + h*(2825.0*k[0][q]/27648.0 +
                18575.0*k[2][q]/48384.0 + 13525.0*k[3][q]/55296.0 +
                277.0*k[4][q]/14336.0 + k[5][q]/4.0);
            error[q] = trial[q] - fourth;
        }
        double sum = 0.0;
#pragma unroll
        for (int q = 0; q < 4; ++q) {
            const double atol = q == 3 ? tolerance*total_speed : tolerance;
            const double scale = atol + tolerance*fmax(fabs(y[q]), fabs(trial[q]));
            const double ratio = error[q]/scale;
            sum = fma(ratio, ratio, sum);
        }
        return sqrt(sum/4.0);
    }
}

template<int Method>
__host__ __device__ __forceinline__ void trace_particle(
    const double *__restrict__ field, const double initial_state[4],
    const Parameters &p, double output[7], unsigned long long counters[3]) {
    double y[4] = {initial_state[0], initial_state[1], initial_state[2],
                   initial_state[3]};
    double first_derivative[4], mu, hmax;
    initialize_particle<Method>(field, p, y, mu, hmax, first_derivative);
    double h = fmax(p.minimum_timestep, 1.0e-3*hmax);
    double hmin_seen = h;
    double hmax_seen = h;
    double t = 0.0;
    double previous_error = 1.0;
    bool first_step = true;
    bool after_reject = false;
    bool lost = hypot(y[0], y[1]) >= 1.0;
    unsigned long long nfev = 1;
    unsigned long long accepted = 0;
    unsigned long long rejected = 0;

    for (int attempt = 0; attempt < kMaximumAttempts && t < p.tmax && !lost;
         ++attempt) {
        h = fmin(h, p.tmax - t);
        double k[kMaximumStages][4];
#pragma unroll
        for (int q = 0; q < 4; ++q) k[0][q] = first_derivative[q];
        double trial[4];
        for (int stage = 1; stage < stage_count<Method>(); ++stage) {
            stage_state<Method>(stage, y, h, k, trial);
            right_hand_side<Method>(field, p, mu, trial, k[stage]);
            ++nfev;
        }
        const double error = finish_attempt<Method>(y, h, k, trial,
                                                     p.tolerance, p.total_speed);
        if constexpr (is_dop853<Method>()) {
            // DOP853's thirteenth RHS evaluation is the endpoint derivative;
            // it is reusable as the first stage after an accepted step.
            right_hand_side<Method>(field, p, mu, trial, k[12]);
            ++nfev;
        }
        const double floor_error = fmax(error, 1.0e-300);
        double factor;
        if constexpr (Method == SIMPLE_CUDA_DORMAND_PRINCE) {
            factor = 0.9*pow(floor_error, -1.0/3.0);
            factor = fmax(0.2, fmin(5.0, factor));
            if (error > 0.5 && error < 1.0) factor = 1.0;
        } else if constexpr (Method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED) {
            if (first_step || after_reject) {
                factor = 0.9*pow(floor_error, -1.0/5.0);
            } else {
                factor = 0.9*pow(floor_error, -0.14)*
                         pow(previous_error, 0.08);
            }
            if (after_reject) factor = fmin(1.0, factor);
            factor = fmax(0.2, fmin(5.0, factor));
        } else if constexpr (is_dop853<Method>()) {
            factor = 0.9*pow(floor_error, -1.0/8.0);
            if (after_reject) factor = fmin(1.0, factor);
            factor = fmax(0.2, fmin(10.0, factor));
        } else if (first_step || after_reject) {
            factor = 0.9*pow(fmax(error, 1.0e-10), -1.0/5.0);
            factor = fmax(0.2, fmin(5.0, factor));
        } else {
            factor = 0.9*pow(fmax(error, 1.0e-10), -0.14)*
                     pow(previous_error, 0.08);
            factor = fmax(0.2, fmin(5.0, factor));
        }
        const double hnew = fmax(p.minimum_timestep, fmin(hmax, h*factor));
        const bool at_minimum = p.minimum_timestep > 0.0 &&
                                h <= p.minimum_timestep;
        if (error <= 1.0 || at_minimum) {
            t += h;
#pragma unroll
            for (int q = 0; q < 4; ++q) y[q] = trial[q];
            hmin_seen = fmin(hmin_seen, h);
            hmax_seen = fmax(hmax_seen, h);
            ++accepted;
            previous_error = fmax(error, 1.0e-10);
            first_step = false;
            after_reject = false;
            lost = hypot(y[0], y[1]) >= 1.0;
            if (!lost && t < p.tmax) {
                if constexpr (is_dopri5<Method>()) {
                    // The seventh DP stage is evaluated at the accepted state,
                    // so it is the first stage of the next attempt (FSAL).
#pragma unroll
                    for (int q = 0; q < 4; ++q)
                        first_derivative[q] = k[6][q];
                } else if constexpr (is_dop853<Method>()) {
#pragma unroll
                    for (int q = 0; q < 4; ++q)
                        first_derivative[q] = k[12][q];
                } else {
                    right_hand_side<Method>(field, p, mu, y, first_derivative);
                    ++nfev;
                }
            }
        } else {
            ++rejected;
            after_reject = true;
        }
        h = hnew;
    }

    output[0] = t;
    output[1] = y[0];
    output[2] = y[1];
    output[3] = y[2];
    output[4] = y[3];
    output[5] = hmin_seen;
    output[6] = hmax_seen;
    counters[0] = nfev;
    counters[1] = accepted;
    counters[2] = rejected;
}

template<int Method>
__global__ void trace_kernel(const double *__restrict__ field,
                             const double *__restrict__ initial_state,
                             Parameters p, int particle_count,
                             double *__restrict__ output,
                             unsigned long long *__restrict__ counters) {
    const int particle = blockIdx.x*blockDim.x + threadIdx.x;
    if (particle >= particle_count) return;
    trace_particle<Method>(field, initial_state + 4*particle, p,
                           output + 7*particle, counters + 3*particle);
}

template<int Method>
void trace_particles_cpu(const double *__restrict__ field,
                         const double *__restrict__ initial_state,
                         const Parameters &p, int particle_count,
                         double *__restrict__ output,
                         unsigned long long *__restrict__ counters) {
#pragma omp parallel for schedule(dynamic, 1)
    for (int particle = 0; particle < particle_count; ++particle)
        trace_particle<Method>(field, initial_state + 4*particle, p,
                               output + 7*particle, counters + 3*particle);
}

using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

int cuda_status(cudaError_t status) {
    return status == cudaSuccess ? 0 : 1000 + static_cast<int>(status);
}

} // namespace

extern "C" int simple_cuda_firm_rk54(
    int method, int particle_count, const double *quad_points,
    size_t quad_point_count, const double ranges[9], const double *initial_stz,
    const double *initial_vparallel, double mass, double charge,
    double total_speed, double psi0, double tmax, double tolerance,
    double minimum_timestep, double *final_stzv,
    unsigned long long *counters, double profile_ms[SIMPLE_CUDA_PROFILE_COUNT]) {
    if (profile_ms) std::fill(profile_ms, profile_ms + SIMPLE_CUDA_PROFILE_COUNT, 0.0);
    if (!supported_method(method)) return 1;
    if (particle_count <= 0 || !quad_points || !ranges || !initial_stz ||
        !initial_vparallel || !final_stzv || !counters || total_speed <= 0.0 ||
        psi0 == 0.0 || tmax < 0.0 || tolerance <= 0.0 ||
        minimum_timestep < 0.0) return 2;

    Parameters p{};
    std::size_t expected = kStoredQuantities*64;
    for (int d = 0; d < 3; ++d) {
        p.minimum[d] = ranges[3*d];
        p.period_or_maximum[d] = ranges[3*d + 1];
        p.point_count[d] = static_cast<int>(ranges[3*d + 2]);
        if (p.point_count[d] < 4 || (p.point_count[d] - 1)%3 != 0 ||
            p.period_or_maximum[d] <= p.minimum[d]) return 3;
        p.step[d] = (p.period_or_maximum[d] - p.minimum[d])/
                    (p.point_count[d] - 1);
        p.cell_count[d] = (p.point_count[d] - 1)/3;
        expected *= static_cast<std::size_t>(p.cell_count[d]);
    }
    if (quad_point_count != expected) return 4;
    p.mass = mass;
    p.charge = charge;
    p.total_speed = total_speed;
    p.psi0 = psi0;
    p.tmax = tmax;
    p.tolerance = tolerance;
    p.minimum_timestep = minimum_timestep;

    const auto total_begin = Clock::now();
    const auto allocate_begin = total_begin;
    std::vector<double> initial_state(4*particle_count);
    for (int particle = 0; particle < particle_count; ++particle) {
        const double s = initial_stz[3*particle];
        const double theta = initial_stz[3*particle + 1];
        initial_state[4*particle] = s*cos(theta);
        initial_state[4*particle + 1] = s*sin(theta);
        initial_state[4*particle + 2] = initial_stz[3*particle + 2];
        initial_state[4*particle + 3] = initial_vparallel[particle];
    }
    double *field_device = nullptr;
    double *state_device = nullptr;
    double *output_device = nullptr;
    unsigned long long *counter_device = nullptr;
    cudaError_t error = cudaMalloc(&field_device, quad_point_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&state_device, 4*particle_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&output_device, 7*particle_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&counter_device, 3*particle_count*sizeof(unsigned long long));
    const auto allocate_end = Clock::now();
    if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_ALLOCATE] =
        milliseconds(allocate_begin, allocate_end);
    if (error != cudaSuccess) goto cleanup;

    {
        const auto upload_begin = Clock::now();
        error = cudaMemcpy(field_device, quad_points,
                           quad_point_count*sizeof(double), cudaMemcpyHostToDevice);
        if (error == cudaSuccess)
            error = cudaMemcpy(state_device, initial_state.data(),
                               4*particle_count*sizeof(double), cudaMemcpyHostToDevice);
        const auto upload_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_UPLOAD] =
            milliseconds(upload_begin, upload_end);
    }
    if (error != cudaSuccess) goto cleanup;

    {
        cudaEvent_t start{}, stop{};
        error = cudaEventCreate(&start);
        if (error == cudaSuccess) error = cudaEventCreate(&stop);
        if (error == cudaSuccess) error = cudaEventRecord(start);
        const int blocks = (particle_count + kThreads - 1)/kThreads;
        if (error == cudaSuccess) {
            if (method == SIMPLE_CUDA_DORMAND_PRINCE) {
                trace_kernel<SIMPLE_CUDA_DORMAND_PRINCE><<<blocks, kThreads>>>(
                    field_device, state_device, p, particle_count,
                    output_device, counter_device);
            } else if (method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED) {
                trace_kernel<SIMPLE_CUDA_DORMAND_PRINCE_TUNED><<<blocks, kThreads>>>(
                    field_device, state_device, p, particle_count,
                    output_device, counter_device);
            } else if (method == SIMPLE_CUDA_DOP853) {
                trace_kernel<SIMPLE_CUDA_DOP853><<<blocks, kThreads>>>(
                    field_device, state_device, p, particle_count,
                    output_device, counter_device);
            } else {
                trace_kernel<SIMPLE_CUDA_CASH_KARP><<<blocks, kThreads>>>(
                    field_device, state_device, p, particle_count,
                    output_device, counter_device);
            }
            error = cudaGetLastError();
        }
        if (error == cudaSuccess) error = cudaEventRecord(stop);
        if (error == cudaSuccess) error = cudaEventSynchronize(stop);
        float kernel_ms = 0.0f;
        if (error == cudaSuccess) error = cudaEventElapsedTime(&kernel_ms, start, stop);
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_KERNEL] = kernel_ms;
        if (start) cudaEventDestroy(start);
        if (stop) cudaEventDestroy(stop);
    }
    if (error != cudaSuccess) goto cleanup;

    {
        const auto download_begin = Clock::now();
        error = cudaMemcpy(final_stzv, output_device,
                           7*particle_count*sizeof(double), cudaMemcpyDeviceToHost);
        if (error == cudaSuccess)
            error = cudaMemcpy(counters, counter_device,
                               3*particle_count*sizeof(unsigned long long),
                               cudaMemcpyDeviceToHost);
        const auto download_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_DOWNLOAD] =
            milliseconds(download_begin, download_end);
    }
    if (error == cudaSuccess) {
        const auto metric_begin = Clock::now();
        for (int particle = 0; particle < particle_count; ++particle) {
            const double x = final_stzv[7*particle + 1];
            const double y = final_stzv[7*particle + 2];
            final_stzv[7*particle + 1] = hypot(x, y);
            final_stzv[7*particle + 2] = atan2(y, x);
        }
        const auto metric_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_METRIC] =
            milliseconds(metric_begin, metric_end);
    }

cleanup:
    cudaFree(counter_device);
    cudaFree(output_device);
    cudaFree(state_device);
    cudaFree(field_device);
    if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_TOTAL] =
        milliseconds(total_begin, Clock::now());
    return cuda_status(error);
}

extern "C" int simple_cuda_firm_rk54_landreman(
    int method, int particle_count, const double *quad_points,
    size_t quad_point_count, const double ranges[9], const double *initial_stz,
    const double *initial_vparallel, double mass, double charge,
    double total_speed, double psi0, double tmax, double tolerance,
    double minimum_timestep, double maxloss, double t_block, double tau,
    double *final_stzv, unsigned long long *counters,
    double profile_ms[SIMPLE_CUDA_PROFILE_COUNT]) {
    if (profile_ms)
        std::fill(profile_ms, profile_ms + SIMPLE_CUDA_PROFILE_COUNT, 0.0);
    if (!supported_method(method)) return 1;
    if (particle_count <= 0 || !quad_points || !ranges || !initial_stz ||
        !initial_vparallel || !final_stzv || !counters || total_speed <= 0.0 ||
        psi0 == 0.0 || tmax <= 0.0 || tolerance <= 0.0 ||
        minimum_timestep < 0.0 || maxloss < 0.0 || maxloss >= 1.0 ||
        t_block <= 0.0 || tau <= 0.0) return 2;

    Parameters p{};
    std::size_t expected = kStoredQuantities*64;
    for (int d = 0; d < 3; ++d) {
        p.minimum[d] = ranges[3*d];
        p.period_or_maximum[d] = ranges[3*d + 1];
        p.point_count[d] = static_cast<int>(ranges[3*d + 2]);
        if (p.point_count[d] < 4 || (p.point_count[d] - 1)%3 != 0 ||
            p.period_or_maximum[d] <= p.minimum[d]) return 3;
        p.step[d] = (p.period_or_maximum[d] - p.minimum[d])/
                    (p.point_count[d] - 1);
        p.cell_count[d] = (p.point_count[d] - 1)/3;
        expected *= static_cast<std::size_t>(p.cell_count[d]);
    }
    if (quad_point_count != expected) return 4;
    p.mass = mass;
    p.charge = charge;
    p.total_speed = total_speed;
    p.psi0 = psi0;
    p.tolerance = tolerance;
    p.minimum_timestep = minimum_timestep;

    const auto total_begin = Clock::now();
    const auto allocate_begin = total_begin;
    std::vector<double> current_state(4*particle_count);
    for (int particle = 0; particle < particle_count; ++particle) {
        const double s = initial_stz[3*particle];
        const double theta = initial_stz[3*particle + 1];
        current_state[4*particle] = s*cos(theta);
        current_state[4*particle + 1] = s*sin(theta);
        current_state[4*particle + 2] = initial_stz[3*particle + 2];
        current_state[4*particle + 3] = initial_vparallel[particle];
    }
    std::vector<double> segment_output(7*particle_count);
    std::vector<unsigned long long> segment_counters(3*particle_count);
    std::vector<double> ordered_state(4*particle_count);
    std::vector<int> particle_order;
    particle_order.reserve(particle_count);
    std::vector<unsigned long long> previous_work(particle_count, 0);
    std::vector<double> loss_time(particle_count, -1.0);
    std::vector<double> timestep_min(particle_count,
                                     std::numeric_limits<double>::infinity());
    std::vector<double> timestep_max(particle_count, 0.0);
    std::vector<unsigned char> lost(particle_count, 0);
    std::fill(counters, counters + 3*particle_count, 0);

    double *field_device = nullptr;
    double *state_device = nullptr;
    double *output_device = nullptr;
    unsigned long long *counter_device = nullptr;
    cudaError_t error = cudaMalloc(&field_device, quad_point_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&state_device, 4*particle_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&output_device, 7*particle_count*sizeof(double));
    if (error == cudaSuccess)
        error = cudaMalloc(&counter_device,
                           3*particle_count*sizeof(unsigned long long));
    const auto allocate_end = Clock::now();
    if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_ALLOCATE] =
        milliseconds(allocate_begin, allocate_end);
    if (error != cudaSuccess) goto segmented_cleanup;

    {
        const auto upload_begin = Clock::now();
        error = cudaMemcpy(field_device, quad_points,
                           quad_point_count*sizeof(double), cudaMemcpyHostToDevice);
        const auto upload_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_UPLOAD] +=
            milliseconds(upload_begin, upload_end);
    }
    if (error != cudaSuccess) goto segmented_cleanup;

    {
        double current_time = 0.0;
        while (current_time < tmax) {
            const double duration = fmin(t_block, tmax - current_time);
            p.tmax = duration;
            particle_order.clear();
            for (int particle = 0; particle < particle_count; ++particle)
                if (!lost[particle]) particle_order.push_back(particle);
            std::stable_sort(
                particle_order.begin(), particle_order.end(),
                [&](int left, int right) {
                    if (previous_work[left] != previous_work[right])
                        return previous_work[left] > previous_work[right];
                    return left > right;
                });
            const int active_count = static_cast<int>(particle_order.size());
            if (active_count == 0) break;
            for (int slot = 0; slot < active_count; ++slot) {
                const int particle = particle_order[slot];
                for (int q = 0; q < 4; ++q)
                    ordered_state[4*slot + q] = current_state[4*particle + q];
            }
            {
                const auto upload_begin = Clock::now();
                error = cudaMemcpy(state_device, ordered_state.data(),
                    4*active_count*sizeof(double), cudaMemcpyHostToDevice);
                const auto upload_end = Clock::now();
                if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_UPLOAD] +=
                    milliseconds(upload_begin, upload_end);
            }
            if (error != cudaSuccess) break;

            const int blocks = (active_count + kThreads - 1)/kThreads;
            cudaEvent_t start{}, stop{};
            error = cudaEventCreate(&start);
            if (error == cudaSuccess) error = cudaEventCreate(&stop);
            if (error == cudaSuccess) error = cudaEventRecord(start);
            if (error == cudaSuccess) {
                if (method == SIMPLE_CUDA_DORMAND_PRINCE) {
                    trace_kernel<SIMPLE_CUDA_DORMAND_PRINCE><<<blocks, kThreads>>>(
                        field_device, state_device, p,
                        active_count, output_device, counter_device);
                } else if (method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED) {
                    trace_kernel<SIMPLE_CUDA_DORMAND_PRINCE_TUNED><<<blocks, kThreads>>>(
                        field_device, state_device, p,
                        active_count, output_device, counter_device);
                } else if (method == SIMPLE_CUDA_DOP853) {
                    trace_kernel<SIMPLE_CUDA_DOP853><<<blocks, kThreads>>>(
                        field_device, state_device, p,
                        active_count, output_device, counter_device);
                } else {
                    trace_kernel<SIMPLE_CUDA_CASH_KARP><<<blocks, kThreads>>>(
                        field_device, state_device, p,
                        active_count, output_device, counter_device);
                }
                error = cudaGetLastError();
            }
            if (error == cudaSuccess) error = cudaEventRecord(stop);
            if (error == cudaSuccess) error = cudaEventSynchronize(stop);
            float kernel_ms = 0.0f;
            if (error == cudaSuccess)
                error = cudaEventElapsedTime(&kernel_ms, start, stop);
            if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_KERNEL] += kernel_ms;
            if (start) cudaEventDestroy(start);
            if (stop) cudaEventDestroy(stop);
            if (error != cudaSuccess) break;

            {
                const auto download_begin = Clock::now();
                error = cudaMemcpy(segment_output.data(), output_device,
                    7*active_count*sizeof(double), cudaMemcpyDeviceToHost);
                if (error == cudaSuccess)
                    error = cudaMemcpy(segment_counters.data(), counter_device,
                        3*active_count*sizeof(unsigned long long),
                        cudaMemcpyDeviceToHost);
                const auto download_end = Clock::now();
                if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_DOWNLOAD] +=
                    milliseconds(download_begin, download_end);
            }
            if (error != cudaSuccess) break;

            const auto metric_begin = Clock::now();
            for (int slot = 0; slot < active_count; ++slot) {
                const int particle = particle_order[slot];
                previous_work[particle] = segment_counters[3*slot];
                timestep_min[particle] = fmin(timestep_min[particle],
                    segment_output[7*slot + 5]);
                timestep_max[particle] = fmax(timestep_max[particle],
                    segment_output[7*slot + 6]);
                for (int c = 0; c < 3; ++c)
                    counters[3*particle + c] += segment_counters[3*slot + c];
                for (int q = 0; q < 4; ++q)
                    current_state[4*particle + q] = segment_output[7*slot + 1 + q];
                if (segment_output[7*slot] < duration - 1.0e-12) {
                    lost[particle] = 1;
                    loss_time[particle] = current_time +
                                          segment_output[7*slot];
                }
            }
            current_time += duration;
            double weighted_losses = 0.0;
            for (int particle = 0; particle < particle_count; ++particle)
                if (lost[particle])
                    weighted_losses += exp(-loss_time[particle]/tau);
            const bool should_stop = weighted_losses/particle_count > maxloss;
            const auto metric_end = Clock::now();
            if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_METRIC] +=
                milliseconds(metric_begin, metric_end);
            if (should_stop) break;
        }

        if (error == cudaSuccess) {
            const auto metric_begin = Clock::now();
            for (int particle = 0; particle < particle_count; ++particle) {
                const double x = current_state[4*particle];
                const double y = current_state[4*particle + 1];
                final_stzv[7*particle] = lost[particle] ?
                    loss_time[particle] : tmax;
                final_stzv[7*particle + 1] = hypot(x, y);
                final_stzv[7*particle + 2] = atan2(y, x);
                final_stzv[7*particle + 3] = current_state[4*particle + 2];
                final_stzv[7*particle + 4] = current_state[4*particle + 3];
                final_stzv[7*particle + 5] = timestep_min[particle];
                final_stzv[7*particle + 6] = timestep_max[particle];
            }
            const auto metric_end = Clock::now();
            if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_METRIC] +=
                milliseconds(metric_begin, metric_end);
        }
    }

segmented_cleanup:
    cudaFree(counter_device);
    cudaFree(output_device);
    cudaFree(state_device);
    cudaFree(field_device);
    if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_TOTAL] =
        milliseconds(total_begin, Clock::now());
    return cuda_status(error);
}

extern "C" int simple_cpu_firm_rk54_landreman(
    int method, int particle_count, const double *quad_points,
    size_t quad_point_count, const double ranges[9], const double *initial_stz,
    const double *initial_vparallel, double mass, double charge,
    double total_speed, double psi0, double tmax, double tolerance,
    double minimum_timestep, double maxloss, double t_block, double tau,
    double *final_stzv, unsigned long long *counters,
    double profile_ms[SIMPLE_CUDA_PROFILE_COUNT]) {
    if (profile_ms)
        std::fill(profile_ms, profile_ms + SIMPLE_CUDA_PROFILE_COUNT, 0.0);
    if (!supported_method(method)) return 1;
    if (particle_count <= 0 || !quad_points || !ranges || !initial_stz ||
        !initial_vparallel || !final_stzv || !counters || total_speed <= 0.0 ||
        psi0 == 0.0 || tmax <= 0.0 || tolerance <= 0.0 ||
        minimum_timestep < 0.0 || maxloss < 0.0 || maxloss >= 1.0 ||
        t_block <= 0.0 || tau <= 0.0) return 2;

    Parameters p{};
    std::size_t expected = kStoredQuantities*64;
    for (int d = 0; d < 3; ++d) {
        p.minimum[d] = ranges[3*d];
        p.period_or_maximum[d] = ranges[3*d + 1];
        p.point_count[d] = static_cast<int>(ranges[3*d + 2]);
        if (p.point_count[d] < 4 || (p.point_count[d] - 1)%3 != 0 ||
            p.period_or_maximum[d] <= p.minimum[d]) return 3;
        p.step[d] = (p.period_or_maximum[d] - p.minimum[d])/
                    (p.point_count[d] - 1);
        p.cell_count[d] = (p.point_count[d] - 1)/3;
        expected *= static_cast<std::size_t>(p.cell_count[d]);
    }
    if (quad_point_count != expected) return 4;
    p.mass = mass;
    p.charge = charge;
    p.total_speed = total_speed;
    p.psi0 = psi0;
    p.tolerance = tolerance;
    p.minimum_timestep = minimum_timestep;

    const auto total_begin = Clock::now();
    const auto allocate_begin = total_begin;
    std::vector<double> current_state(4*particle_count);
    for (int particle = 0; particle < particle_count; ++particle) {
        const double s = initial_stz[3*particle];
        const double theta = initial_stz[3*particle + 1];
        current_state[4*particle] = s*cos(theta);
        current_state[4*particle + 1] = s*sin(theta);
        current_state[4*particle + 2] = initial_stz[3*particle + 2];
        current_state[4*particle + 3] = initial_vparallel[particle];
    }
    std::vector<double> segment_output(7*particle_count);
    std::vector<unsigned long long> segment_counters(3*particle_count);
    std::vector<double> loss_time(particle_count, -1.0);
    std::vector<double> timestep_min(particle_count,
                                     std::numeric_limits<double>::infinity());
    std::vector<double> timestep_max(particle_count, 0.0);
    std::vector<unsigned char> lost(particle_count, 0);
    std::fill(counters, counters + 3*particle_count, 0);
    const auto allocate_end = Clock::now();
    if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_ALLOCATE] =
        milliseconds(allocate_begin, allocate_end);

    double current_time = 0.0;
    while (current_time < tmax) {
        const double duration = fmin(t_block, tmax - current_time);
        p.tmax = duration;
        const auto trace_begin = Clock::now();
        if (method == SIMPLE_CUDA_DORMAND_PRINCE) {
            trace_particles_cpu<SIMPLE_CUDA_DORMAND_PRINCE>(
                quad_points, current_state.data(), p, particle_count,
                segment_output.data(), segment_counters.data());
        } else if (method == SIMPLE_CUDA_DORMAND_PRINCE_TUNED) {
            trace_particles_cpu<SIMPLE_CUDA_DORMAND_PRINCE_TUNED>(
                quad_points, current_state.data(), p, particle_count,
                segment_output.data(), segment_counters.data());
        } else if (method == SIMPLE_CUDA_DOP853) {
            trace_particles_cpu<SIMPLE_CUDA_DOP853>(
                quad_points, current_state.data(), p, particle_count,
                segment_output.data(), segment_counters.data());
        } else {
            trace_particles_cpu<SIMPLE_CUDA_CASH_KARP>(
                quad_points, current_state.data(), p, particle_count,
                segment_output.data(), segment_counters.data());
        }
        const auto trace_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_KERNEL] +=
            milliseconds(trace_begin, trace_end);

        const auto metric_begin = Clock::now();
        for (int particle = 0; particle < particle_count; ++particle) {
            timestep_min[particle] = fmin(timestep_min[particle],
                segment_output[7*particle + 5]);
            timestep_max[particle] = fmax(timestep_max[particle],
                segment_output[7*particle + 6]);
            for (int c = 0; c < 3; ++c)
                counters[3*particle + c] += segment_counters[3*particle + c];
            if (lost[particle]) continue;
            current_state[4*particle] = segment_output[7*particle + 1];
            current_state[4*particle + 1] = segment_output[7*particle + 2];
            current_state[4*particle + 2] = segment_output[7*particle + 3];
            current_state[4*particle + 3] = segment_output[7*particle + 4];
            if (segment_output[7*particle] < duration - 1.0e-12) {
                lost[particle] = 1;
                loss_time[particle] = current_time +
                                      segment_output[7*particle];
            }
        }
        current_time += duration;
        double weighted_losses = 0.0;
        for (int particle = 0; particle < particle_count; ++particle)
            if (lost[particle])
                weighted_losses += exp(-loss_time[particle]/tau);
        const bool should_stop = weighted_losses/particle_count > maxloss;
        const auto metric_end = Clock::now();
        if (profile_ms) profile_ms[SIMPLE_CUDA_PROFILE_METRIC] +=
            milliseconds(metric_begin, metric_end);
        if (should_stop) break;
    }

    const auto metric_begin = Clock::now();
    for (int particle = 0; particle < particle_count; ++particle) {
        const double x = current_state[4*particle];
        const double y = current_state[4*particle + 1];
        final_stzv[7*particle] = lost[particle] ? loss_time[particle] : tmax;
        final_stzv[7*particle + 1] = hypot(x, y);
        final_stzv[7*particle + 2] = atan2(y, x);
        final_stzv[7*particle + 3] = current_state[4*particle + 2];
        final_stzv[7*particle + 4] = current_state[4*particle + 3];
        final_stzv[7*particle + 5] = timestep_min[particle];
        final_stzv[7*particle + 6] = timestep_max[particle];
    }
    const auto metric_end = Clock::now();
    if (profile_ms) {
        profile_ms[SIMPLE_CUDA_PROFILE_METRIC] +=
            milliseconds(metric_begin, metric_end);
        profile_ms[SIMPLE_CUDA_PROFILE_TOTAL] =
            milliseconds(total_begin, Clock::now());
    }
    return 0;
}

extern "C" const char *simple_cuda_error_string(int status) {
    switch (status) {
    case 0: return "success";
    case 1: return "unsupported integrator";
    case 2: return "invalid argument";
    case 3: return "invalid FIRM3D metagrid range";
    case 4: return "FIRM3D metagrid size mismatch";
    default:
        if (status >= 1000)
            return cudaGetErrorString(static_cast<cudaError_t>(status - 1000));
        return "unknown SIMPLE CUDA error";
    }
}
