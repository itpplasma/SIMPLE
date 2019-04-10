#include <stdio.h>

#include "../matlab/neo_orb.h"

int main(int argc, char **argv) {
  int ns_s = 3;
  int ns_tp = 3;
  int multharm = 5;
  int integmode = 1;
  __neo_orb_MOD_init_field(&ns_s, &ns_tp, &multharm, &integmode);
}
