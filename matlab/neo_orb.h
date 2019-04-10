#include <stddef.h>

void __neo_orb_MOD_init_field(int *ans_s, int *ans_tp, int *amultharm,
  int *aintegmode);
void __neo_orb_MOD_init_params(int *Z_charge, int *m_mass, double *E_kin,
  double *adtau, double *adtaumax, double *arelerr);
void __neo_orb_MOD_timestep(double *s, double *th, double *ph, double *lam,
  int *ierr);
void magfie_vmec_(double *x, double *bmod, double *sqrtg, double *bder,
  double *hcovar, double *hctrvr, double *hcurl);
void magfie_can_(double *x, double *bmod, double *sqrtg, double *bder,
  double *hcovar, double *hctrvr, double *hcurl);
