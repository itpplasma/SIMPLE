import numpy as np
import matplotlib.pyplot as plt
import gvec


state = gvec.load_state(
    "build/test/tests/gvec_vmec_params.ini",
    "build/test/tests/test_vmec_gvec_State_0000_00000000.dat",
)
ev = state.evaluate("mod_B", rho=20, theta=20, zeta=0.0)

ev.mod_B.plot()
