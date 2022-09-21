"""
Created: Fri Aug  9 15:50:40 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""
#%%

import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from bench_interface import *

nskip = 20

runtimes = []
evalss = []
Herr = []
jparerr = []
runtimes_bench = {}
runtimes_single = {}
evals_bench = {}
evals_single = {}
Herr_bench = {}
avgdist_single = {}
Herr_single = {}
jparerr_bench = {}

runtime = None
evals = None
#%%
bench.multi = 0
bench.quasi = 0
bench.tok = 1
bench.integ_mode = 1
bench.npoiper2 = 500000
bench.nlag = 0
bench.nplagr_invar = 4
bench.ncut = 10000
bench.nt = 500000
bench.rtol = 1e-13
bench.init_bench()


#%% Getting reference
    
bench.nt = 500000
bench.npoiper2 = 500000
bench.ncut = 0
bench.multi = 0
integ_mode = 3
bench.do_bench()
refdata = np.loadtxt('/tmp/out.txt')
xref = 1.0 + refdata[:,0]*np.cos(refdata[:,1])
yref = refdata[:,0]*np.sin(refdata[:,1])
zref = np.stack((xref, yref)).copy(order='F')

bench.nt = 10000
res = do_singlebench(4, [3], nskip, zref)
print(res)
data = np.loadtxt('/tmp/out.txt')[::nskip,:]
x = 1.0 + data[:,0]*np.cos(data[:,1])
y = data[:,0]*np.sin(data[:,1])
plt.figure()
plt.plot(x, y, ',')
plt.plot(xref, yref, ',')
#bench.cleanup_bench()  
#import scipy.interpolate as spi
#tck, u = spi.splprep(refdata[:,:2].T, s=0, k=3)
#
#def refsqdist(t, z):
#    splres = spi.splev(t, tck)
#    return (z[0]-splres[0])**2 + (z[1]-splres[1])**2


#%% Benchmarking Hamiltonian and distance
bench.nt = 80000
bench.npoiper2 = 8

bench.quasi = 1
bench.integ_mode = 0
bench.ncut = 0
tols = np.array([1e-6, 1e-8, 1e-10, 1e-12])
avgdist_single[0] = []
runtimes_single[0] = []
evals_single[0] = []
Herr_single[0] = []
for tol in tols:
    bench.rtol = tol
    bench.do_bench()
    data = np.loadtxt('/tmp/out.txt')[::nskip,:]
    x = 1.0 + data[:,0]*np.cos(data[:,1])
    y = data[:,0]*np.sin(data[:,1])
    avgdist_single[0].append(avgdist(np.stack((x,y)).copy(order='F'), zref))
    Hmean = np.sqrt(np.mean(data[1:,4]**2))
    Herr_single[0].append(np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2)))
    runtimes_single[0].append(bench.endtime - bench.starttime)
    evals_single[0].append(diag.icounter)
evals_single[0] = np.array(evals_single[0])
Herr_single[0] = np.array(Herr_single[0])
runtimes_single[0] = np.array(runtimes_single[0])
avgdist_single[0] = np.array(avgdist_single[0])
#%%
data = np.loadtxt('/tmp/out.txt')[::nskip,:]
x = 1.0 + data[:,0]*np.cos(data[:,1])
y = data[:,0]*np.sin(data[:,1])
plt.figure()
plt.plot(x, y, ',')
plt.plot(xref, yref, ',')
#%%
bench.nt = 8000
bench.quasi = 0
bench.rtol = 1e-13

npoipers = np.array([3, 8, 16, 32, 64, 128, 256], int)
integ_mode = 1
bench.nlag = 1
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:], nskip, zref)
integ_mode = 2
bench.nlag = 0
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:], nskip, zref)
integ_mode = 3
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode], 
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers, nskip, zref)
integ_mode = 4
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers, nskip, zref)

bench.multi = 1
integ_mode = 21
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode], 
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:], nskip, zref)
    
#%%
plt.figure()
plt.loglog(np.array(runtimes_single[1])*npoipers[1:]/bench.nt, avgdist_single[1], 'k-')
plt.loglog(np.array(runtimes_single[0])/10000, avgdist_single[0], '-', color='tab:gray')
plt.loglog(np.array(runtimes_single[3])*npoipers/bench.nt, avgdist_single[3], 'k:')
plt.loglog(np.array(runtimes_single[4])*npoipers/bench.nt, avgdist_single[4], 'k--')
plt.loglog(np.array(runtimes_single[21])*npoipers[1:]/bench.nt, avgdist_single[21], 'k-.')
#plt.loglog(np.array(runtimes_single[5])*npoipers/bench.nt, avgdist_single[5], 'k-.')
#plt.loglog(np.array(runtimes_single[15])*npoipers[1:]/bench.nt, avgdist_single[15], 'r-.')
plt.xlim([3e-6, 1e-3])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('CPU time per bounce period / s')
plt.ylabel('$\delta x_{\mathrm{pol}}$')

exportpng('fig_tok_dist1')

plt.figure()
plt.loglog(np.array(evals_single[1])*npoipers[1:]/bench.nt, avgdist_single[1], 'k-')
plt.loglog(np.array(evals_single[0])/10000, avgdist_single[0], '-', color='tab:gray')
plt.loglog(np.array(evals_single[3])*npoipers/bench.nt, avgdist_single[3], 'k:')
plt.loglog(np.array(evals_single[4])*npoipers/bench.nt, avgdist_single[4], 'k--')
plt.loglog(np.array(evals_single[21])*npoipers[1:]/bench.nt, avgdist_single[21], 'k-.')
#plt.loglog(np.array(evals_single[5])*npoipers/bench.nt, avgdist_single[5], 'k-.')
#plt.loglog(np.array(evals_single[15])*npoipers[1:]/bench.nt, avgdist_single[15], 'r-.')
plt.xlim([1e1, 1.1e3])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('evaluations $N_{\mathrm{f}}$ per bounce period')
plt.ylabel('$\delta x_{\mathrm{pol}}$')

exportpng('fig_tok_dist2')

#%%
plt.figure()
plt.loglog(np.array(runtimes_single[1])*npoipers[1:]/bench.nt, Herr_single[1], 'k-')
plt.loglog(np.array(runtimes_single[0])/10000, Herr_single[0], '-', color='tab:gray')
plt.loglog(np.array(runtimes_single[3])*npoipers/bench.nt, Herr_single[3], 'k:')
plt.loglog(np.array(runtimes_single[4])*npoipers/bench.nt, Herr_single[4], 'k--')
plt.loglog(np.array(runtimes_single[21])*npoipers[1:]/bench.nt, Herr_single[21], 'k-.')
#plt.loglog(np.array(runtimes_single[5])*npoipers/bench.nt, Herr_single[5], 'k-.')
#plt.loglog(np.array(runtimes_single[15])*npoipers[1:]/bench.nt, Herr_single[15], 'r-.')
#plt.xlim([4e-3, 2e0])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('CPU time / s')
plt.ylabel('$\delta H$')

plt.figure()
plt.loglog(np.array(evals_single[1])*npoipers[1:]/bench.nt, Herr_single[1], 'k-')
plt.loglog(np.array(evals_single[0])/10000, Herr_single[0], '-', color='tab:gray')
plt.loglog(np.array(evals_single[3])*npoipers/bench.nt, Herr_single[3], 'k:')
plt.loglog(np.array(evals_single[4])*npoipers/bench.nt, Herr_single[4], 'k--')
plt.loglog(np.array(evals_single[21])*npoipers[1:]/bench.nt, Herr_single[21], 'k-.')
#plt.loglog(np.array(evals_single[5])*npoipers/bench.nt, Herr_single[5], 'k-.')
#plt.loglog(np.array(evals_single[15])*npoipers[1:]/bench.nt, Herr_single[15], 'r-.')
#plt.xlim([1e4, 1e7])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')


#%%
