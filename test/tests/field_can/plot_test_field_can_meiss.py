import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('lam_chi.out')
R = data[:,0]
PH = data[:,1]
TH = data[:,2]
lam = data[:,3]
chi = data[:,4]

r = np.unique(R)
ph = np.unique(PH)
th = np.unique(TH)


plt.figure()
plt.title('lam')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], lam[condi])

plt.figure()
plt.title('chi')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], chi[condi])

#%%
data = np.loadtxt('covar_components.out')
R = data[:,0]
PH = data[:,1]
TH = data[:,2]
Ar = data[:,3]
Ap = data[:,4]
At = data[:,5]
hr = data[:,6]
hp = data[:,7]
ht = data[:,8]

r = np.unique(R)
ph = np.unique(PH)
th = np.unique(TH)

plt.figure()
plt.title('Ar')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], Ar[condi])

plt.figure()
plt.title('Ap')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], Ap[condi])

plt.figure()
plt.title('hr')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], hr[condi])

plt.figure()
plt.title('hp')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], hp[condi])

dz1 = -hr/hp
dz2 = Ar + Ap*dz1

plt.figure()
plt.title('dz1')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], dz1[condi])

plt.figure()
plt.title('dz2')
for phk in ph:
    for thk in th:
        condi = np.logical_and(np.isclose(TH, thk), np.isclose(PH, phk))
        plt.plot(R[condi], dz2[condi])

plt.show()
