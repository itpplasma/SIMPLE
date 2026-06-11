#include <netcdf.h>
#include <stddef.h>
#include <stdio.h>

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
        n_rho = 3,
        n_theta = 4,
        n_zeta = 5,
        n_theta_field = n_theta + 1,
        n_zeta_field = n_zeta + 1
    };

    const char *filename = argc > 1 ? argv[1] : "legacy_boozer_chartmap.nc";
    double rho[n_rho] = {0.001, 0.5005, 1.0};
    double s[n_rho] = {0.000001, 0.5000005, 1.0};
    double theta[n_theta] = {0.0, 1.5707963267948966, 3.141592653589793,
                             4.71238898038469};
    double zeta[n_zeta] = {0.0, 0.3141592653589793, 0.6283185307179586,
                           0.9424777960769379, 1.2566370614359172};
    double a_phi[n_rho] = {0.0, 0.1, 0.2};
    double b_theta[n_rho] = {1.0, 1.0, 1.0};
    double b_phi[n_rho] = {2.0, 2.0, 2.0};
    double xyz[(size_t)n_zeta * n_theta * n_rho];
    double bmod[(size_t)n_zeta_field * n_theta_field * n_rho];
    int ncid;
    int dim_rho, dim_s, dim_theta, dim_zeta, dim_theta_field, dim_zeta_field;
    int dims_geom[3], dims_field[3];
    int var_rho, var_s, var_theta, var_zeta, var_x, var_y, var_z;
    int var_aphi, var_btheta, var_bphi, var_bmod, var_nfp;

    for (size_t i = 0; i < sizeof(xyz) / sizeof(xyz[0]); ++i) {
        xyz[i] = 1.0;
    }
    for (size_t i = 0; i < sizeof(bmod) / sizeof(bmod[0]); ++i) {
        bmod[i] = 1.0;
    }

    CHECK_NC(nc_create(filename, NC_NETCDF4, &ncid), "create file");
    CHECK_NC(nc_def_dim(ncid, "rho", n_rho, &dim_rho), "def dim rho");
    CHECK_NC(nc_def_dim(ncid, "s", n_rho, &dim_s), "def dim s");
    CHECK_NC(nc_def_dim(ncid, "theta", n_theta, &dim_theta), "def dim theta");
    CHECK_NC(nc_def_dim(ncid, "zeta", n_zeta, &dim_zeta), "def dim zeta");
    CHECK_NC(nc_def_dim(ncid, "theta_field", n_theta_field, &dim_theta_field),
             "def dim theta_field");
    CHECK_NC(nc_def_dim(ncid, "zeta_field", n_zeta_field, &dim_zeta_field),
             "def dim zeta_field");

    dims_geom[0] = dim_zeta;
    dims_geom[1] = dim_theta;
    dims_geom[2] = dim_rho;
    dims_field[0] = dim_zeta_field;
    dims_field[1] = dim_theta_field;
    dims_field[2] = dim_rho;

    CHECK_NC(nc_def_var(ncid, "rho", NC_DOUBLE, 1, &dim_rho, &var_rho), "def rho");
    CHECK_NC(nc_def_var(ncid, "s", NC_DOUBLE, 1, &dim_s, &var_s), "def s");
    CHECK_NC(nc_def_var(ncid, "theta", NC_DOUBLE, 1, &dim_theta, &var_theta),
             "def theta");
    CHECK_NC(nc_def_var(ncid, "zeta", NC_DOUBLE, 1, &dim_zeta, &var_zeta),
             "def zeta");
    CHECK_NC(nc_def_var(ncid, "x", NC_DOUBLE, 3, dims_geom, &var_x), "def x");
    CHECK_NC(nc_def_var(ncid, "y", NC_DOUBLE, 3, dims_geom, &var_y), "def y");
    CHECK_NC(nc_def_var(ncid, "z", NC_DOUBLE, 3, dims_geom, &var_z), "def z");
    CHECK_NC(nc_def_var(ncid, "A_phi", NC_DOUBLE, 1, &dim_s, &var_aphi),
             "def A_phi");
    CHECK_NC(nc_def_var(ncid, "B_theta", NC_DOUBLE, 1, &dim_rho, &var_btheta),
             "def B_theta");
    CHECK_NC(nc_def_var(ncid, "B_phi", NC_DOUBLE, 1, &dim_rho, &var_bphi),
             "def B_phi");
    CHECK_NC(nc_def_var(ncid, "Bmod", NC_DOUBLE, 3, dims_field, &var_bmod),
             "def Bmod");
    CHECK_NC(nc_def_var(ncid, "num_field_periods", NC_INT, 0, NULL, &var_nfp),
             "def num_field_periods");
    CHECK_NC(nc_put_att_double(ncid, NC_GLOBAL, "torflux", NC_DOUBLE, 1,
                               (const double[]){1.0}),
             "put torflux");
    CHECK_NC(nc_put_att_text(ncid, var_aphi, "radial_abscissa", 1, "s"),
             "put A_phi radial_abscissa");
    CHECK_NC(nc_enddef(ncid), "enddef");

    CHECK_NC(nc_put_var_double(ncid, var_rho, rho), "put rho");
    CHECK_NC(nc_put_var_double(ncid, var_s, s), "put s");
    CHECK_NC(nc_put_var_double(ncid, var_theta, theta), "put theta");
    CHECK_NC(nc_put_var_double(ncid, var_zeta, zeta), "put zeta");
    CHECK_NC(nc_put_var_double(ncid, var_x, xyz), "put x");
    CHECK_NC(nc_put_var_double(ncid, var_y, xyz), "put y");
    CHECK_NC(nc_put_var_double(ncid, var_z, xyz), "put z");
    CHECK_NC(nc_put_var_double(ncid, var_aphi, a_phi), "put A_phi");
    CHECK_NC(nc_put_var_double(ncid, var_btheta, b_theta), "put B_theta");
    CHECK_NC(nc_put_var_double(ncid, var_bphi, b_phi), "put B_phi");
    CHECK_NC(nc_put_var_double(ncid, var_bmod, bmod), "put Bmod");
    CHECK_NC(nc_put_var_int(ncid, var_nfp, (const int[]){5}), "put num_field_periods");
    CHECK_NC(nc_close(ncid), "close");
    return 0;
}
