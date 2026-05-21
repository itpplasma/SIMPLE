# %%
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root
from common import (dt_cpp_var, r0, th0,ph0,vpar0,timesteps,autosave,folderforPlots,mu,qe,m,c,)
from field_correct_test import field


f = field()

nt = timesteps


metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}

def midpoint_quantities(x, xold, dt):
    xmid = 0.5*(x + xold)
    vmid = (x - xold)/dt
    f.evaluate(xmid[0], xmid[1], xmid[2])
    Amid = np.array([0, f.co_Ath, f.co_Aph])
    g = metric(xmid)
    dLdx = (
        m/2*(g['d_11']*vmid[0]**2 + g['d_22']*vmid[1]**2 + g['d_33']*vmid[2]**2)
        + qe/c*(f.co_dAth*vmid[1] + f.co_dAph*vmid[2])
        - mu*f.dB
    )
    dLdxdot = m*g['_ii']*vmid + qe/c*Amid
    return dLdx, dLdxdot, xmid, vmid


def F(x, xold, dLdxold, dLdxdotold, dt):
    dLdx, dLdxdot, _, _ = midpoint_quantities(x, xold, dt)
    return (dLdx + dLdxold)*dt/2 - (dLdxdot - dLdxdotold)

''' --- IGNORE ---
def F(x, xold, dLdxold, dLdxdotold, dt):
    global p, dpdt
    ret = np.zeros(3)

    xmid = (x + xold)/2
    vmid = (x - xold)/dt

    f.evaluate(xmid[0], xmid[1], xmid[2])
    Amid = np.array([0, f.co_Ath, f.co_Aph])
    g = metric(xmid)

    dpdt =  m/2*(g['d_11']*vmid[0]**2 + g['d_22']*vmid[1]**2 + g['d_33']*vmid[2]**2) + qe/c*(f.co_dAth*vmid[1] + f.co_dAph*vmid[2]) - mu*f.dB
    p = m*g['_ii']*vmid + qe/c*Amid

    ret = (dpdt + dLdxold)*dt/2 - (p - dLdxdotold)
    return ret
'''


Z = []
Vel = []
Vparvec = []
Vpervec = []
Vpar = []
Vper = []
V = []
P = []
Energy = []
Pphi = []
Mu = []
from time import time
for i in range(len(dt_cpp_var)):
    dt = dt_cpp_var[i]
    #nt = int(nt / dt)

    z = np.zeros([3, nt + 1])
    z[:, 0] = [r0, th0, ph0]

    vel = np.zeros([3, nt + 1])

    f.evaluate(r0, th0, ph0)
    g = metric(z[:,0])

    A = np.array([0, f.co_Ath, f.co_Aph])
    dLdxdotold = m*g['_ii']*vel[:,0] + qe/c*A
    dLdxold = m/2*(g['d_11']*vel[0,0]**2 + g['d_22']*vel[1,0]**2 + g['d_33']*vel[2,0]**2) + qe/c*(f.co_dAth*vel[1,0] + f.co_dAph*vel[2,0]) - mu*f.dB

    P_con = np.zeros([3, nt+1])
    P_con[:,0] = g['^ii']*dLdxdotold.copy()

    v_per_vec = np.zeros([3, nt+1])
    v_par_vec = np.zeros([3, nt+1])
    v_par_vec[:,0] = vel[:,0]
    vper = np.zeros(nt+1)
    vpar = np.zeros(nt+1)
    v = np.zeros(nt+1)

    H = np.zeros(nt+1)
    pphi = np.zeros(nt+1)
    MU = np.zeros(nt+1)
    H[0] = mu * f.B
    pphi[0] = m*g['_ii'][2]*vel[2,0] + qe/c*f.co_Aph
    MU[0] = mu


    tic = time()
    for kt in range(nt):
        p = m*g['_ii']*vel[:,0] + qe/c*A

        sol = root(F, z[:,kt], method='hybr',tol=1e-12,args=(z[:,kt], dLdxold, dLdxdotold, dt))
        z[:,kt+1] = sol.x

        dLdxold, dLdxdotold, xmid, vmid = midpoint_quantities(z[:,kt+1], z[:,kt], dt)
        f.evaluate(z[0,kt+1], z[1,kt+1], z[2,kt+1])
        g = metric(z[:,kt+1])
        A = np.array([0, f.co_Ath, f.co_Aph])
        vel[:,kt+1] = 1/m * g['^ii']*((p + dt*dLdxold) - qe/c*A)



        P_con[:,kt+1] = g['^ii']*p.copy()

        v_par_vec[:,kt+1] = [0,
                            g['^ii'][1]*f.co_hth*(vel[1,kt+1]*f.co_hth + vel[2,kt+1]*f.co_hph),
                            g['^ii'][2]*f.co_hph*(vel[2,kt+1]*f.co_hph + vel[1,kt+1]*f.co_hth)]
        v_per_vec[:,kt+1] = vel[:,kt+1] - v_par_vec[:,kt+1]

        scalarvper =np.sqrt((v_per_vec[0,kt+1]**2 * g['_ii'][0] + v_per_vec[1,kt+1]**2 * g['_ii'][1] + v_per_vec[2,kt+1]**2 * g['_ii'][2]))
        scalarvpar =np.sqrt((v_par_vec[0,kt+1]**2 * g['_ii'][0] + v_par_vec[1,kt+1]**2 * g['_ii'][1] + v_par_vec[2,kt+1]**2 * g['_ii'][2]))


        vper[kt+1] = scalarvper #np.sqrt(g['^ii'][0]*mu*2*f.B)
        vpar[kt+1] = scalarvpar #np.sqrt(g['_ii'][0]*vel[0,kt+1]**2 + g['_ii'][1]*vel[1,kt+1]**2 + g['_ii'][2]*vel[2,kt+1]**2) - vper[kt]
        v[kt+1] = np.sqrt((vel[0,kt+1]**2 * g['_ii'][0] + vel[1,kt+1]**2 * g['_ii'][1] + vel[2,kt+1]**2 * g['_ii'][2]))

        H[kt+1] = 1/(2*m)*(g['^ii'][0]*p[0]**2 + g['^ii'][1]*(p[1] - qe/c*f.co_Ath)**2 + g['^ii'][2]*(p[2] - qe/c*f.co_Aph)**2) + mu*f.B
        pphi[kt+1] = m*g['_ii'][2]*vel[2,kt+1]+ qe/c*f.co_Aph
        MU[kt+1] = m*vper[kt+1]**2/(2*f.B)

    Z.append(z)
    Vel.append(vel)
    Vparvec.append(v_par_vec)
    Vpervec.append(v_per_vec)
    Vpar.append(vpar)
    Vper.append(vper)
    V.append(v)
    P.append(P_con)
    Energy.append(H)
    Pphi.append(pphi)
    Mu.append(MU)


q_cpp_var = Z
v_cpp_var = Vel
p_cpp_var = P
Vtot_cpp_var = V
Vpar_cpp_var = Vpar
Vper_cpp_var = Vper
vpar_cpp_var_vec = Vparvec
vper_cpp_var_vec = Vpervec
En_cpp_var = Energy
Pphi_cpp_var = Pphi
Mu_cpp_var = Mu

k = nt#//100
fig = plt.figure(figsize=(18,6))

for j in range(len(dt_cpp_var)):
    t = np.linspace(0, nt*dt_cpp_var[j], nt)
    ax = fig.add_subplot(1, len(dt_cpp_var), j+1)
    ax.plot(t[:k], Vpar[j][:k], linewidth=0.5, c='r', label=r'$v_\parallel$')
    ax.plot(t[:k], Vper[j][:k], linewidth=0.5, c='g', label=r'$v_\perp$')
    ax.plot(t[:k], V[j][:k], linewidth=0.5, c='b', label=r'$v$')
    ax.legend()
    ax.set_title(r'$\Delta t_{\rm cpp\_var} = ' + str(dt_cpp_var[j]) + '$')

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'vel_cpp_var.png'), dpi=300)
else:
    plt.show()




fig = plt.figure(figsize=(18,6))

for j in range(len(dt_cpp_var)):
    ax = fig.add_subplot(1, len(dt_cpp_var), j+1)
    ax.scatter(Z[j][0,:]*np.cos(Z[j][1,:]), Z[j][0,:]*np.sin(Z[j][1,:]), s=0.01, color='g')
    ax.set_title(r'$\Delta t_{\rm cpp\_var} = ' + str(dt_cpp_var[j]) + '$')

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'banana_cpp_var.png'), dpi=300)
plt.show()
