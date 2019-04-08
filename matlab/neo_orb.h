#include <stddef.h>

void __neo_orb_MOD_init_field(int *ans_s, int *ans_tp, int *amultharm);
void __neo_orb_MOD_init_params(int *Z_charge, int *m_mass, double *E_kin,
  double *adtau, double *adtaumax, int *aintegmode, double *arelerr);
void __neo_orb_MOD_timestep(double *s, double *th, double *ph, double *lam,
  int *ierr);
