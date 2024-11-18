
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import f, r0, th0, ph0, pph0, timesteps, get_val, get_der
from plotting import plot_orbit, plot_cost_function

dt, nt = timesteps(steps_per_bounce = 8, nbounce = 100)
nlag = 1 # order of Lagrange extrapolation

qe = 1.0; m = 1.0; c = 1.0 # particle charge, mass and speed of light
mu = 1e-5 # magnetic moment

z = np.zeros([3,nt+1])
z[:,0] = [r0, th0, ph0]
[H, pth, vpar] = get_val(np.array([r0,th0,ph0,pph0]))


calc_fields = lambda vpar: {
    'A2': m * vpar * f.hth + qe * f.Ath / c,
    'dA2': np.array([
        m * vpar * f.dhth[0] + qe * f.dAth[0] / c,
        m * vpar * f.dhth[1] + qe * f.dAth[1] / c,
        m * vpar * f.dhth[2] + qe * f.dAth[2] / c
    ]),
    'A3': m * vpar * f.hph + qe * f.Aph / c,
    'dA3': np.array([
        m * vpar * f.dhph[0] + qe * f.dAph[0] / c,
        m * vpar * f.dhph[1] + qe * f.dAph[1] / c,
        m * vpar * f.dhph[2] + qe * f.dAph[2] / c
    ]),
    'dB': mu * f.dB,
    'dPhie': f.dPhie
}


def F(x, zold, eq2, eq3):
    f.evaluate(x[0], x[1], zold[2])
    cov = calc_fields(x[2])

    eq1 = -(cov['dA2'][0] *(x[1] - zold[1]) - dt *(cov['dB'][0] + cov['dPhie'][0])) / (cov['dA3'][0]) 

    ret = np.zeros(3)
    ret[0] = eq2 - (m*x[2]*f.hth + qe*f.Ath/c)
    ret[1] = eq3 - (m*x[2]*f.hph + qe*f.Aph/c)
    ret[2] = m*f.hth *(x[1] - zold[1]) + m*f.hph*eq1 - dt*x[2]

    return ret



from time import time
tic = time()
for kt in range(nt):
    vparold = vpar

    # Initialize via Lagrange extrapolation
    if(kt>=nlag):
        x0 = np.zeros(2)
        extrapr0 = lagrange(np.arange(-nlag, 1), z[0,kt-nlag:kt+1])
        extrapr1 = lagrange(np.arange(-nlag, 1), z[1,kt-nlag:kt+1])
        x0[0] = extrapr0(1.0)
        x0[1] = extrapr1(1.0)
    else:
        x0 = np.array([z[0,kt], z[1,kt]])


    f.evaluate(z[0,kt], z[1,kt], z[2,kt])  
    cov = calc_fields(vparold)
    # M_ij * delta_j = b_i
    M = np.array([[cov['dA2'][0], cov['dA3'][0]],[m*f.hth,  m*f.hph]])
    b = dt * np.array([(cov['dB'][0] + cov['dPhie'][0]), (vparold)])
    Delta = np.linalg.solve(M, b)

    eq2 = cov['dA2'][1]*Delta[0] + cov['dA3'][1]*Delta[1] + cov['A2'] - dt*(cov['dB'][1] + cov['dPhie'][1])
    eq3 = cov['dA2'][2]*Delta[0] + cov['dA3'][2]*Delta[1] + cov['A3'] - dt*(cov['dB'][2] + cov['dPhie'][2])
    
    # Implicit substep for A+
    X0 = np.append(x0 , vparold)
    sol = root(F, X0, method='hybr',tol=1e-12,args=(z[:,kt], eq2, eq3))
    r, th, vpar = sol.x
    
    f.evaluate(r, th, z[2,kt+1])
    cov = calc_fields(vpar)
    # M_ij * delta_j = b_i
    M = np.array([[cov['dA2'][0], cov['dA3'][0]],[m*f.hth,  m*f.hph]])
    b = dt * np.array([(cov['dB'][0] + cov['dPhie'][0]), (vpar)])
    Delta = np.linalg.solve(M, b)

    
    z[0,kt+1] = r
    z[1,kt+1] = z[1,kt] + Delta[0]
    z[2,kt+1] = 0



    
plot_orbit(z)
plt.show()

