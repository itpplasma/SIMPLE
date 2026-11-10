/* Counter-based standard normal numbers for collision steps.
 *
 * The number for (seed, marker id, step, component) depends on nothing else,
 * so SIMPLE and ASCOT5 (itpplasma/ascot5 src/simulate/mccc/collision_noise.h,
 * which this file copies; keep both identical) apply identical collision kicks
 * to the same marker when they share the seed and the fixed step.
 *
 * Mixing is splitmix64 (Steele, Lea, Flood 2014); the normal deviate is
 * Box-Muller from two 53-bit uniforms.
 */
#include <math.h>
#include <stdint.h>

static inline uint64_t collision_noise_mix(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline double collision_noise_normal(uint64_t seed, uint64_t id,
                                            uint64_t step, uint64_t component) {
    uint64_t key = collision_noise_mix(
        seed ^ collision_noise_mix(
            id ^ collision_noise_mix(step ^ collision_noise_mix(component))));
    uint64_t a = collision_noise_mix(key);
    uint64_t b = collision_noise_mix(key ^ 0xd1b54a32d192ed03ULL);
    double u1 = ((double)(a >> 11) + 0.5) * 0x1.0p-53;
    double u2 = (double)(b >> 11) * 0x1.0p-53;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

double simple_collision_noise_normal(int64_t seed, int64_t id, int64_t step,
                                     int64_t component) {
    return collision_noise_normal((uint64_t)seed, (uint64_t)id, (uint64_t)step,
                                  (uint64_t)component);
}
