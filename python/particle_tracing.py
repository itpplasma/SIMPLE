
from pysimple import simple, simple_main, params, orbit_symplectic, params
from pysimple import get_can_sub as coord
from pysimple import field_can_mod

import numpy as np
import matplotlib.pyplot as plt
#%%
tracy = simple.Tracer()

params.read_config("simple.in")
simple_main.init_field(tracy, "wout.nc", 3, 3, 3, 1)
params.params_init()

field_can_mod.eval_field(tracy.f, tracy.si.z[0], tracy.si.z[1], tracy.si.z[2])
