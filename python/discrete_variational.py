# %%
import os
import matplotlib.pyplot as plt
import numpy as np
from common import (dt_gc_var, f,r0,th0,ph0,pph0,timesteps,autosave,folderforPlots,get_val,mu,qe,m,c,)
from field_correct_test import field



f = field()
nt = timesteps

print('nt =', nt)



# metric is only needed for the calculation of velocities for the plots!
metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}


calc_fields = lambda vpar: {
    "A2": m * vpar * f.co_hth + qe * f.co_Ath / c,
    "dA2": m * vpar * f.co_dhth + qe * f.co_dAth / c,
    "d2A2": m * vpar * f.co_d2hth + qe * f.co_d2Ath / c,
    "A3": m * vpar * f.co_hph + qe * f.co_Aph / c,
    "dA3": m * vpar * f.co_dhph + qe * f.co_dAph / c,
    "d2A3": m * vpar * f.co_d2hph + qe * f.co_d2Aph / c,
    "dB": mu * f.dB,
    "d2B": mu * f.d2B,
    "dPhie": f.dPhie,
    "d2Phie": f.d2Phie,
}


def Newton(x, zold, eq2, eq3, dt):

    eps = 1
    while eps >= 1e-11:
        f.evaluate(x[0], x[1], zold[2])
        cov = calc_fields(x[2])

        eq1 = -(cov["dA2"][0] * (x[1] - zold[1]) - dt * (cov["dB"][0] + cov["dPhie"][0])) / (cov["dA3"][0])
        deq1 = np.array([
                - (cov["d2A2"][0] * (x[1] - zold[1]) - dt * (cov["d2B"][0] + cov["d2Phie"][0])) / (cov["dA3"][0])
                + (cov["dA2"][0] * (x[1] - zold[1]) - dt * (cov["dB"][0] + cov["dPhie"][0])) / (cov["dA3"][0]) ** 2
                * cov["d2A3"][0],
                - (cov["d2A2"][3] * (x[1] - zold[1]) + cov["dA2"][0] - dt * (cov["d2B"][3] + cov["d2Phie"][3])) / (cov["dA3"][0])
                + (cov["dA2"][0] * (x[1] - zold[1]) - dt * (cov["dB"][0] + cov["dPhie"][0])) / (cov["dA3"][0]) ** 2
                * cov["d2A3"][3],
                - (m * f.co_dhth[0] * (x[1] - zold[1])) / (cov["dA3"][0])
                + (cov["dA2"][0] * (x[1] - zold[1]) - dt * (cov["dB"][0] + cov["dPhie"][0])) / (cov["dA3"][0]) ** 2
                * (m * f.co_dhph[0]),
            ])

        F = np.zeros(3)
        F[0] = eq2 - cov["A2"]
        F[1] = eq3 - cov["A3"]
        F[2] = m * f.co_hth * (x[1] - zold[1]) + m * f.co_hph * eq1 - dt * x[2]

        Jacobi = np.array([
                [-cov["dA2"][0], -cov["dA2"][1], -m * f.co_hth],
                [-cov["dA3"][0], -cov["dA3"][1], -m * f.co_hph],
                [m * (f.co_dhth[0] * (x[1] - zold[1]) + f.co_dhph[0] * eq1 + f.co_hph * deq1[0]),
                 m * (f.co_hth + f.co_dhph[1] * eq1 + f.co_hph * deq1[1]),
                 m * f.co_hph * deq1[2] - dt],
        ])

        delta = np.linalg.solve(Jacobi, -F)

        xnew = x + delta
        x = xnew

        eps = abs(F[2])

    return x



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
for i in range(len(dt_gc_var)):
    dt = dt_gc_var[i]
    #nt = int(nt / dt)

    f.evaluate(r0, th0, ph0)

    z = np.zeros([3, nt + 1])
    z[:, 0] = [r0, th0, ph0]
    [H, pth, vpar] = get_val(np.array([r0, th0, ph0, pph0]))

    vel = np.zeros([3, nt + 1])
    P_con = np.zeros([3, nt+1])
    v_per_vec = np.zeros([3, nt+1])
    v_par_vec = np.zeros([3, nt+1])
    v_par_vec[:,0] = vel[:,0]

    vper = np.zeros(nt+1)
    vparr = np.zeros(nt+1)
    v = np.zeros(nt+1)
    H = np.zeros(nt+1)
    pphi = np.zeros(nt+1)
    MU = np.zeros(nt+1)

    H[0] =  mu*f.B
    pphi[0] = qe/c*f.co_Aph
    MU[0] = mu
    v[0] = np.sum(vel[:,0])
    vper[0] = 0
    vparr[0] = vpar

    tic = time()
    for kt in range(nt):
        vparold = vpar

        f.evaluate(z[0, kt], z[1, kt], z[2, kt])
        g = metric(z[:,kt])
        cov = calc_fields(vparold)

        p = np.array([0, cov["A2"], cov["A3"]])
        P_con[:,kt] = g['^ii'] * p



        # M_ij * delta_j = b_i
        M = np.array([[cov["dA2"][0], cov["dA3"][0]], [m * f.co_hth, m * f.co_hph]])
        b = dt * np.array([(cov["dB"][0] + cov["dPhie"][0]), (vparold)])
        Delta = np.linalg.solve(M, b)

        eq2 = (
            cov["dA2"][1] * Delta[0]
            + cov["dA3"][1] * Delta[1]
            + cov["A2"]
            - dt * (cov["dB"][1] + cov["dPhie"][1])
        )
        eq3 = (
            cov["dA2"][2] * Delta[0]
            + cov["dA3"][2] * Delta[1]
            + cov["A3"]
            - dt * (cov["dB"][2] + cov["dPhie"][2])
        )

        # Implicit substep with Newton
        x0 = np.array([z[0, kt], z[1, kt], vparold])
        x = Newton(x0, z[:, kt], eq2, eq3, dt)

        r = x[0]
        th = x[1]
        vpar = x[2]

        f.evaluate(r, th, z[2, kt + 1])
        cov = calc_fields(vpar)
        # M_ij * delta_j = b_i
        M = np.array([[cov["dA2"][0], cov["dA3"][0]], [m * f.co_hth, m * f.co_hph]])
        b = dt * np.array([(cov["dB"][0] + cov["dPhie"][0]), (vpar)])
        Delta = np.linalg.solve(M, b)

        z[0, kt + 1] = r
        z[1, kt + 1] = z[1, kt] + Delta[0]
        z[2, kt + 1] = z[2, kt] + Delta[1]


        f.evaluate(z[0,kt+1], z[1,kt+1], z[2,kt+1])
        g = metric(z[:,kt+1])
        vel[:,kt+1] = 1/m * g['^ii']*(np.array([0, cov["A2"], cov["A3"]]) - qe/c*np.array([0, f.co_Ath, f.co_Aph]))

        v_par_vec[:,kt+1] = [0,
                            g['^ii'][1]*f.co_hth*(vel[1,kt+1]*f.co_hth + vel[2,kt+1]*f.co_hph),
                            g['^ii'][2]*f.co_hph*(vel[2,kt+1]*f.co_hph + vel[1,kt+1]*f.co_hth)]
        v_per_vec[:,kt+1] = vel[:,kt+1] - v_par_vec[:,kt+1]

        scalarvper =np.sqrt((v_per_vec[0,kt+1]**2 * g['_ii'][0] + v_per_vec[1,kt+1]**2 * g['_ii'][1] + v_per_vec[2,kt+1]**2 * g['_ii'][2]))
        scalarvpar =np.sqrt((v_par_vec[0,kt+1]**2 * g['_ii'][0] + v_par_vec[1,kt+1]**2 * g['_ii'][1] + v_par_vec[2,kt+1]**2 * g['_ii'][2]))

        vper[kt+1] = scalarvper #np.sqrt(g['^ii'][0]*mu*2*f.B)
        vparr[kt+1] = scalarvpar #np.sqrt(np.sum(vel[:,kt+1]**2)) - vper[kt]
        v[kt+1] = np.sqrt((vel[0,kt+1]**2 * g['_ii'][0] + vel[1,kt+1]**2 * g['_ii'][1] + vel[2,kt+1]**2 * g['_ii'][2]))

        H[kt+1] = 1/(2*m)*(g['^ii'][0]*p[0]**2 + g['^ii'][1]*(p[1] - qe/c*f.co_Ath)**2 + g['^ii'][2]*(p[2] - qe/c*f.co_Aph)**2) + mu*f.B
        pphi[kt+1] = m*g['_ii'][2]*vel[2,kt+1]+ qe/c*f.co_Aph
        MU[kt+1] = m*vper[kt+1]**2/(2*f.B)




    Z.append(z)
    Vel.append(vel)
    Vparvec.append(v_par_vec)
    Vpervec.append(v_per_vec)
    Vpar.append(vparr)
    Vper.append(vper)
    V.append(v)
    P.append(P_con)
    Energy.append(H)
    Pphi.append(pphi)
    Mu.append(MU)


q_gc = Z
v_gc = Vel
p_gc = P
Vtot_gc = V
Vpar_gc = Vpar
Vper_gc = Vper
vpar_gc_vec = Vparvec
vper_gc_vec = Vpervec
En_gc = Energy
Pphi_gc = Pphi
Mu_gc = Mu

k = nt#//100
fig = plt.figure(figsize=(18,6))

for j in range(len(dt_gc_var)):
    t = np.linspace(0, nt*dt_gc_var[j], nt)
    ax = fig.add_subplot(1, len(dt_gc_var), j+1)
    ax.plot(t[:k], Vpar[j][:k], linewidth=0.5, c='r', label=r'$v_\parallel$')
    ax.plot(t[:k], Vper[j][:k], linewidth=0.5, c='g', label=r'$v_\perp$')
    ax.plot(t[:k], V[j][:k], linewidth=0.5, c='b', label=r'$v$')
    ax.legend()
    ax.set_title(r'$\Delta t_{\rm gc\_var} = ' + str(dt_gc_var[j]) + '$')

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'vel_gc_var.png'), dpi=300)
else:
    plt.show()



fig = plt.figure(figsize=(18,6))

for j in range(len(dt_gc_var)):
    ax = fig.add_subplot(1, len(dt_gc_var), j+1)
    ax.scatter(Z[j][0,:]*np.cos(Z[j][1,:]), Z[j][0,:]*np.sin(Z[j][1,:]), s=0.01, color='purple')
    ax.set_title(r'$\Delta t_{\rm gc\_var} = ' + str(dt_gc_var[j]) + '$')

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'banana_gc_var.png'), dpi=300)
plt.show()
