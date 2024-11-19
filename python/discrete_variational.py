#%%
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
    'dA2': m * vpar * f.dhth + qe * f.dAth / c,
    'd2A2': m * vpar * f.d2hth + qe * f.d2Ath / c,
    'A3': m * vpar * f.hph + qe * f.Aph / c,
    'dA3': m * vpar * f.dhph + qe * f.dAph / c,
    'd2A3': m * vpar * f.d2hph + qe * f.d2Aph /c,

    'dB': mu * f.dB,
    'd2B': mu * f.d2B,
    'dPhie': f.dPhie,
    'd2Phie': f.d2Phie
}


def Newton(x, zold, eq2, eq3):

    eps = 1
    while eps >= 1e-11:
        f.evaluate(x[0], x[1], zold[2])
        cov = calc_fields(x[2])

        eq1 = -(cov['dA2'][0] *(x[1] - zold[1]) - dt *(cov['dB'][0] + cov['dPhie'][0])) / (cov['dA3'][0]) 
        deq1 = np.array([-(cov['d2A2'][0] *(x[1]-zold[1]) - dt*(cov['d2B'][0] + cov['d2Phie'][0])) / (cov['dA3'][0]) + (cov['dA2'][0] *(x[1] - zold[1]) - dt *(cov['dB'][0] + cov['dPhie'][0])) /(cov['dA3'][0])**2 * cov['d2A3'][0],
                     -(cov['d2A2'][3] *(x[1]-zold[1]) + cov['dA2'][0] - dt*(cov['d2B'][3] + cov['d2Phie'][3])) / (cov['dA3'][0]) + (cov['dA2'][0] *(x[1] - zold[1]) - dt *(cov['dB'][0] + cov['dPhie'][0])) /(cov['dA3'][0])**2 * cov['d2A3'][3],
                     -(m*f.dhth[0] *(x[1]- zold[1]))/(cov['dA3'][0]) + (cov['dA2'][0] *(x[1] - zold[1]) - dt *(cov['dB'][0] + cov['dPhie'][0])) /(cov['dA3'][0])**2 * (m * f.dhph[0])])

        F = np.zeros(3)
        F[0] = eq2 - cov['A2']
        F[1] = eq3 - cov['A3']
        F[2] = m*f.hth *(x[1] - zold[1]) + m*f.hph*eq1 - dt*x[2]

        Jacobi = np.array([[-cov['dA2'][0], -cov['dA2'][1], -m*f.hth], 
                      [-cov['dA3'][0], -cov['dA3'][1], -m*f.hph],
                      [m*(f.dhth[0]*(x[1] - zold[1]) + f.dhph[0]*eq1 + f.hph*deq1[0]), m*(f.hth + f.dhph[1]*eq1 + f.hph*deq1[1]), m*f.hph*deq1[2] - dt]])

        delta = np.linalg.solve(Jacobi, -F)

        xnew = x + delta
        x = xnew

        eps = abs(F[2])



    return x



from time import time
tic = time()
for kt in range(nt):
    vparold = vpar


    f.evaluate(z[0,kt], z[1,kt], z[2,kt])  
    cov = calc_fields(vparold)
    # M_ij * delta_j = b_i
    M = np.array([[cov['dA2'][0], cov['dA3'][0]],[m*f.hth,  m*f.hph]])
    b = dt * np.array([(cov['dB'][0] + cov['dPhie'][0]), (vparold)])
    Delta = np.linalg.solve(M, b)

    eq2 = cov['dA2'][1]*Delta[0] + cov['dA3'][1]*Delta[1] + cov['A2'] - dt*(cov['dB'][1] + cov['dPhie'][1])
    eq3 = cov['dA2'][2]*Delta[0] + cov['dA3'][2]*Delta[1] + cov['A3'] - dt*(cov['dB'][2] + cov['dPhie'][2])
    
    # Implicit substep with Newton
    x0 = np.array([z[0,kt], z[1,kt], vparold])
    x = Newton(x0, z[:,kt], eq2, eq3)

    r = x[0]
    th = x[1]
    vpar = x[2]
    
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

# %%
