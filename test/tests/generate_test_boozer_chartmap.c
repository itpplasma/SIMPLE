#include <math.h>
#include <netcdf.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define CHECK_NC(call, context)                                                   \
    do {                                                                          \
        int status_ = (call);                                                     \
        if (status_ != NC_NOERR) {                                                \
            fprintf(stderr, "NetCDF error in %s: %s\n", context,                  \
                    nc_strerror(status_));                                        \
            return 1;                                                             \
        }                                                                         \
    } while (0)

int main(int argc, char **argv)
{
    enum
    {
        nfp = 2,
        n_rho = 17,
        n_theta = 33,
        n_zeta = 33
    };

    const double b0 = 2.0e4;
    const double r0 = 150.0;
    const double a = 50.0;
    const double iota = 0.5;
    const double twopi = 2.0 * acos(-1.0);
    const double epsilon = a / r0;
    const double torflux = 0.5 * b0 * a * a;
    const char *filename = argc > 1 ? argv[1] : "test_boozer_chartmap.nc";

    double rho[n_rho];
    double theta[n_theta];
    double zeta[n_zeta];
    double a_phi[n_rho];
    double b_theta[n_rho];
    double b_phi[n_rho];
    double *x = NULL;
    double *y = NULL;
    double *z = NULL;
    double *bmod = NULL;
    int ncid;
    int dim_rho, dim_theta, dim_zeta;
    int dims_3d[3];
    int var_rho, var_theta, var_zeta, var_x, var_y, var_z, var_aphi, var_btheta;
    int var_bphi, var_bmod, var_nfp;

    x = malloc((size_t)n_rho * n_theta * n_zeta * sizeof(*x));
    y = malloc((size_t)n_rho * n_theta * n_zeta * sizeof(*y));
    z = malloc((size_t)n_rho * n_theta * n_zeta * sizeof(*z));
    bmod = malloc((size_t)n_rho * n_theta * n_zeta * sizeof(*bmod));
    if (x == NULL || y == NULL || z == NULL || bmod == NULL) {
        fprintf(stderr, "Allocation failure\n");
        free(x);
        free(y);
        free(z);
        free(bmod);
        return 1;
    }

    for (int ir = 0; ir < n_rho; ++ir) {
        rho[ir] = (double)ir / (double)(n_rho - 1);
        a_phi[ir] = -torflux * iota * rho[ir] * rho[ir];
        b_theta[ir] = b0 * epsilon * iota;
        b_phi[ir] = b0 * r0;
    }
    rho[0] = 1.0e-6;
    a_phi[0] = -torflux * iota * rho[0] * rho[0];

    for (int it = 0; it < n_theta; ++it) {
        theta[it] = twopi * (double)it / (double)n_theta;
    }

    for (int iz = 0; iz < n_zeta; ++iz) {
        zeta[iz] = (twopi / (double)nfp) * (double)iz / (double)n_zeta;
    }

    for (int iz = 0; iz < n_zeta; ++iz) {
        const double zeta_val = zeta[iz];
        for (int it = 0; it < n_theta; ++it) {
            const double theta_val = theta[it];
            const double cos_theta = cos(theta_val);
            const double sin_theta = sin(theta_val);
            for (int ir = 0; ir < n_rho; ++ir) {
                const double rho_val = rho[ir];
                const double r_minor = a * rho_val;
                const double radius = r0 + r_minor * cos_theta;
                const size_t idx = ((size_t)iz * n_theta + (size_t)it) * n_rho +
                                   (size_t)ir;
                x[idx] = radius * cos(zeta_val);
                y[idx] = radius * sin(zeta_val);
                z[idx] = r_minor * sin_theta;
                bmod[idx] = b0 / (1.0 + rho_val * epsilon * cos_theta);
            }
        }
    }

    CHECK_NC(nc_create(filename, NC_NETCDF4, &ncid), "create file");
    CHECK_NC(nc_def_dim(ncid, "rho", n_rho, &dim_rho), "def dim rho");
    CHECK_NC(nc_def_dim(ncid, "theta", n_theta, &dim_theta), "def dim theta");
    CHECK_NC(nc_def_dim(ncid, "zeta", n_zeta, &dim_zeta), "def dim zeta");

    dims_3d[0] = dim_zeta;
    dims_3d[1] = dim_theta;
    dims_3d[2] = dim_rho;

    CHECK_NC(nc_def_var(ncid, "rho", NC_DOUBLE, 1, &dim_rho, &var_rho), "def rho");
    CHECK_NC(nc_def_var(ncid, "theta", NC_DOUBLE, 1, &dim_theta, &var_theta), "def theta");
    CHECK_NC(nc_def_var(ncid, "zeta", NC_DOUBLE, 1, &dim_zeta, &var_zeta), "def zeta");
    CHECK_NC(nc_def_var(ncid, "x", NC_DOUBLE, 3, dims_3d, &var_x), "def x");
    CHECK_NC(nc_def_var(ncid, "y", NC_DOUBLE, 3, dims_3d, &var_y), "def y");
    CHECK_NC(nc_def_var(ncid, "z", NC_DOUBLE, 3, dims_3d, &var_z), "def z");
    CHECK_NC(nc_def_var(ncid, "A_phi", NC_DOUBLE, 1, &dim_rho, &var_aphi), "def A_phi");
    CHECK_NC(nc_def_var(ncid, "B_theta", NC_DOUBLE, 1, &dim_rho, &var_btheta), "def B_theta");
    CHECK_NC(nc_def_var(ncid, "B_phi", NC_DOUBLE, 1, &dim_rho, &var_bphi), "def B_phi");
    CHECK_NC(nc_def_var(ncid, "Bmod", NC_DOUBLE, 3, dims_3d, &var_bmod), "def Bmod");
    CHECK_NC(nc_def_var(ncid, "num_field_periods", NC_INT, 0, NULL, &var_nfp),
             "def num_field_periods");

    CHECK_NC(nc_put_att_text(ncid, var_x, "units", 2, "cm"), "put x units");
    CHECK_NC(nc_put_att_text(ncid, var_y, "units", 2, "cm"), "put y units");
    CHECK_NC(nc_put_att_text(ncid, var_z, "units", 2, "cm"), "put z units");
    CHECK_NC(nc_put_att_text(ncid, NC_GLOBAL, "rho_convention", 7, "rho_tor"),
             "put rho_convention");
    CHECK_NC(nc_put_att_text(ncid, NC_GLOBAL, "zeta_convention", 6, "boozer"),
             "put zeta_convention");
    CHECK_NC(nc_put_att_double(ncid, NC_GLOBAL, "rho_lcfs", NC_DOUBLE, 1,
                               (const double[]){1.0}),
             "put rho_lcfs");
    CHECK_NC(nc_put_att_int(ncid, NC_GLOBAL, "boozer_field", NC_INT, 1,
                            (const int[]){1}),
             "put boozer_field");
    CHECK_NC(nc_put_att_double(ncid, NC_GLOBAL, "torflux", NC_DOUBLE, 1, &torflux),
             "put torflux");

    CHECK_NC(nc_enddef(ncid), "enddef");

    CHECK_NC(nc_put_var_double(ncid, var_rho, rho), "put rho");
    CHECK_NC(nc_put_var_double(ncid, var_theta, theta), "put theta");
    CHECK_NC(nc_put_var_double(ncid, var_zeta, zeta), "put zeta");
    CHECK_NC(nc_put_var_double(ncid, var_x, x), "put x");
    CHECK_NC(nc_put_var_double(ncid, var_y, y), "put y");
    CHECK_NC(nc_put_var_double(ncid, var_z, z), "put z");
    CHECK_NC(nc_put_var_double(ncid, var_aphi, a_phi), "put A_phi");
    CHECK_NC(nc_put_var_double(ncid, var_btheta, b_theta), "put B_theta");
    CHECK_NC(nc_put_var_double(ncid, var_bphi, b_phi), "put B_phi");
    CHECK_NC(nc_put_var_double(ncid, var_bmod, bmod), "put Bmod");
    CHECK_NC(nc_put_var_int(ncid, var_nfp, (const int[]){nfp}), "put num_field_periods");
    CHECK_NC(nc_close(ncid), "close");

    printf("Generated %s\n", filename);
    printf("  torflux=%14.6e, B0=%8.1f, R0=%8.1f\n", torflux, b0, r0);
    printf("  nfp=%d, nrho=%d, ntheta=%d, nzeta=%d\n", nfp, n_rho, n_theta, n_zeta);

    free(x);
    free(y);
    free(z);
    free(bmod);
    return 0;
}
