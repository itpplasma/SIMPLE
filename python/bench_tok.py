"""
Created: Fri Aug  9 15:50:40 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

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

#%%

def do_singlebench(integ_mode, npoipers):
    runtimes = []
    evalss = []
    avgdists = []
    Herr = []
    
    bench.ncut = 0
    bench.integ_mode = integ_mode
    for npoiper in npoipers:
        bench.npoiper2 = npoiper
        bench.do_bench()
        data = np.loadtxt('/tmp/out.txt')[::nskip,:]
        x = 1.0 + data[:,0]*np.cos(data[:,1])
        y = data[:,0]*np.sin(data[:,1])
        avgdists.append(avgdist(np.stack((x,y)).copy(order='F')))
        Hmean = np.sqrt(np.mean(data[1:,4]**2))
        Herr.append(np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2)))
        runtimes.append(bench.endtime - bench.starttime)
        evalss.append(diag.icounter)
        
    return np.array(runtimes), np.array(evalss), np.array(avgdists), np.array(Herr)
    
    

def do_run(integ_mode, npoiper2, nplagr_invar):
    global data, runtimes, evalss, Herr, jparerr, jparmean, Hmean, runtime, evals
    
    bench.integ_mode = integ_mode
    bench.npoiper2 = npoiper2
    bench.nplagr_invar = nplagr_invar
    
    bench.do_bench()
    
    runtime = bench.endtime - bench.starttime
    evals = diag.icounter
    
    data = np.loadtxt('/tmp/out.txt')
    Hmean = np.sqrt(np.mean(data[1:,4]**2))
    jparmean = np.sqrt(np.mean(data[1:,5]**2))
    
    runtimes.append(runtime)
    evalss.append(evals)
    Herr.append(np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2)))
    jparerr.append(np.sqrt(np.mean((data[1:,5]/jparmean-1.0)**2)))

def do_plot():
    global data, runtimes, evalss, Herr, jparerr, jparmean, Hmean, runtime, evals
    #plt.figure()
    #plt.plot(data[:,0]*np.cos(data[:,1]), data[:,0]*np.sin(data[:,1]), ',')
    
    plt.figure()
    plt.plot(data[1:,4]/Hmean, '-')
    
    #plt.figure()
    #plt.plot(data[1:,5]/jparmean, '-')
    
def do_cutbench(integ_mode, npoipers, nplagr_invar):
    global data, runtimes, evalss, Herr, jparerr, jparmean, Hmean, runtime, evals
    runtimes = []
    evalss = []
    Herr = []
    jparerr = []
    nrun = len(npoipers)
    print('')
    for krun in range(nrun):
        do_run(integ_mode, npoipers[krun], nplagr_invar[krun])
        #if(krun == 1): do_plot()
        print('{} {:.2e} {} {:.2e} {:.2e}'.format(npoipers[krun],
            runtime, evals, Herr[-1], jparerr[-1]))
    
    return np.array(runtimes), np.array(evalss), np.array(Herr), np.array(jparerr)

#def avgdist(z):
#    nx = z.shape[0]
#    krange = range(nx)    
#    sqdists = np.empty(nx)
#    tlast = 0.0
#    for k in krange:
#        res = spo.minimize(refsqdist, tlast, args = z[k,:], tol=1e-12)
#        #res = spo.minimize_scalar(refsqdist, args = z[k,:], tol=1e-12)
#        sqdists[k] = refsqdist(res.x, z[k,:])
#        #print(res.x, sqdists[k])
#        tlast = res.x + 0.01
#        if tlast > 1: tlast = 0.01
#        
#    #print(np.sqrt(sqdists))
#    return np.sqrt(np.mean(sqdists))

def avgdist(z):
    nz = z.shape[1]
    sqdists = np.empty(nz)
    bench.minsqdist(z, zref, sqdists)
    #print(sqdists)
    return np.sqrt(np.mean(sqdists))

#%% Getting reference
    
bench.nt = 500
bench.npoiper2 = 500
bench.ncut = 0
bench.multi = 0
integ_mode = 3
bench.do_bench()
refdata = np.loadtxt('/tmp/out.txt')
xref = 1.0 + refdata[:,0]*np.cos(refdata[:,1])
yref = refdata[:,0]*np.sin(refdata[:,1])
zref = np.stack((xref, yref)).copy(order='F')

bench.nt = 10000
res = do_singlebench(4, [3])
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
bench.nt = 8000
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
    avgdist_single[0].append(avgdist(np.stack((x,y)).copy(order='F')))
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
bench.quasi = 0
bench.rtol = 1e-13

npoipers = np.array([3, 8, 16, 32, 64, 128, 256], int)
integ_mode = 1
bench.nlag = 1
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:])
integ_mode = 2
bench.nlag = 0
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:])
integ_mode = 3
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode], 
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers)
integ_mode = 4
(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode],
 Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers)
#integ_mode = 5
#(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode], 
# Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers)
#integ_mode = 15
#(runtimes_single[integ_mode], evals_single[integ_mode], avgdist_single[integ_mode], 
# Herr_single[integ_mode]) = do_singlebench(integ_mode, npoipers[1:])
    
#%%
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(np.array(runtimes_single[1])*npoipers[1:]/bench.nt, avgdist_single[1], 'k-')
plt.loglog(np.array(runtimes_single[0])/10000, avgdist_single[0], 'k--')
plt.loglog(np.array(runtimes_single[3])*npoipers/bench.nt, avgdist_single[3], 'k:')
plt.loglog(np.array(runtimes_single[4])*npoipers/bench.nt, avgdist_single[4], 'k-.')
#plt.loglog(np.array(runtimes_single[5])*npoipers/bench.nt, avgdist_single[5], 'k-.')
#plt.loglog(np.array(runtimes_single[15])*npoipers[1:]/bench.nt, avgdist_single[15], 'r-.')
plt.xlim([3e-6, 1e-3])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('CPU time / s')
plt.ylabel('$\delta x$')

plt.subplot(1,2,2)
plt.loglog(np.array(evals_single[1])*npoipers[1:]/bench.nt, avgdist_single[1], 'k-')
plt.loglog(np.array(evals_single[0])/10000, avgdist_single[0], 'k--')
plt.loglog(np.array(evals_single[3])*npoipers/bench.nt, avgdist_single[3], 'k:')
plt.loglog(np.array(evals_single[4])*npoipers/bench.nt, avgdist_single[4], 'k-.')
#plt.loglog(np.array(evals_single[5])*npoipers/bench.nt, avgdist_single[5], 'k-.')
#plt.loglog(np.array(evals_single[15])*npoipers[1:]/bench.nt, avgdist_single[15], 'r-.')
plt.xlim([1e1, 1.1e3])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('field evaluations')
plt.ylabel('$\delta x$')
plt.tight_layout()

#%%
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(np.array(runtimes_single[1])*npoipers[1:]/bench.nt, Herr_single[1], 'k-')
plt.loglog(np.array(runtimes_single[2])*npoipers[1:]/bench.nt, Herr_single[2], 'k--')
plt.loglog(np.array(runtimes_single[3])*npoipers/bench.nt, Herr_single[3], 'k:')
plt.loglog(np.array(runtimes_single[4])*npoipers/bench.nt, Herr_single[4], 'k--')
#plt.loglog(np.array(runtimes_single[5])*npoipers/bench.nt, Herr_single[5], 'k-.')
#plt.loglog(np.array(runtimes_single[15])*npoipers[1:]/bench.nt, Herr_single[15], 'r-.')
#plt.xlim([4e-3, 2e0])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('runtime / s')
plt.ylabel('$\delta H$')

plt.subplot(1,2,2)
plt.loglog(np.array(evals_single[1])*npoipers[1:]/bench.nt, Herr_single[1], 'k-')
plt.loglog(np.array(evals_single[2])*npoipers[1:]/bench.nt, Herr_single[2], 'k--')
plt.loglog(np.array(evals_single[3])*npoipers/bench.nt, Herr_single[3], 'k:')
plt.loglog(np.array(evals_single[4])*npoipers/bench.nt, Herr_single[4], 'k-.')
#plt.loglog(np.array(evals_single[5])*npoipers/bench.nt, Herr_single[5], 'k-.')
#plt.loglog(np.array(evals_single[15])*npoipers[1:]/bench.nt, Herr_single[15], 'r-.')
#plt.xlim([1e4, 1e7])
#plt.ylim([1e-14, 1e-1])
plt.ylim([1e-6, 5e-2])
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')
plt.tight_layout()
    
    
#%% Benchmarking parallel invariant
ncut = 10000
bench.ncut = ncut

# RK45
bench.quasi = 1
bench.npoiper = 16
bench_multi = 0
tols = np.array([5e-7, 1e-8, 1e-10, 1e-12])
print('')
for tol in tols:
    bench.rtol = tol
    do_run(0, 16, 4)
    print('{} {:.2e} {} {:.2e} {:.2e}'.format(4,
        runtime, evals, Herr[-1], jparerr[-1]))

runtimes_bench[0] = np.array(runtimes)
evals_bench[0] = np.array(evalss)
Herr_bench[0] = np.array(Herr)
jparerr_bench[0] = np.array(jparerr)
#%%
bench.quasi = 0
bench.rtol = 1e-13

bench.multi = 1
# Stellarator
#npoipers = [20, 24, 28, 32, 48, 64, 128, 256, 512]

# Tokamak
bench.nt = 10000

npoipers = [32, 64, 128, 256, 512]
nplagr_invar = 6*np.ones(len(npoipers), int)
integ_mode = 22
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

npoipers = [8, 16, 32, 64, 128, 256]
nplagr_invar = np.concatenate([[2,4], 6*np.ones(len(npoipers)-2, int)])

integ_mode = 21
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

bench.multi = 0

integ_mode = 1
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)
    
#integ_mode = 2
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 3
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 4
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 15
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

# Stellarator
#npoipers = [20, 24, 28, 32, 64, 128]
# Tokamak
npoipers = [5, 8, 16, 32, 64, 128, 256]
nplagr_invar = np.concatenate([[2,2,4], 6*np.ones(len(npoipers)-2, int)])
integ_mode = 5
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

#integ_mode = 6
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)
    

#%%
# Tokamak
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[1]/ncut, jparerr_bench[1], 'k-')
plt.loglog(runtimes_bench[3]/ncut, jparerr_bench[3], 'k:')
plt.loglog(runtimes_bench[4]/ncut, jparerr_bench[4], 'k--')
plt.loglog(runtimes_bench[21]/ncut, jparerr_bench[21], 'k-.')
#plt.xlim([4e-3, 2e0])
plt.ylim([1e-13, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta J_{\parallel} $')
plt.subplot(1,2,2)
plt.loglog(evals_bench[1]/ncut, jparerr_bench[1], 'k-')
plt.loglog(evals_bench[3]/ncut, jparerr_bench[3], 'k:')
plt.loglog(evals_bench[4]/ncut, jparerr_bench[4], 'k--')
plt.loglog(evals_bench[21]/ncut, jparerr_bench[21], 'k-.')
#plt.xlim([1e4, 1e7])
plt.ylim([1e-13, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta J_{\parallel} $')
plt.tight_layout()
#%%

plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[0], Herr_bench[0], 'r-')
#%%
plt.loglog(runtimes_bench[1], Herr_bench[1], 'k-')
plt.loglog(runtimes_bench[3], Herr_bench[3], 'k:')
plt.loglog(runtimes_bench[4], Herr_bench[4], 'k--')
plt.loglog(runtimes_bench[21], Herr_bench[21], 'k-.')
plt.xlim([4e-3, 2e0])
plt.ylim([1e-14, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta H$')

plt.subplot(1,2,2)
plt.loglog(evals_bench[1], Herr_bench[1], 'k-')
plt.loglog(evals_bench[3], Herr_bench[3], 'k:')
plt.loglog(evals_bench[4], Herr_bench[4], 'k--')
plt.loglog(evals_bench[21], Herr_bench[21], 'k-.')
plt.xlim([1e4, 1e7])
plt.ylim([1e-14, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')
plt.tight_layout()
    
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[15], jparerr_bench[15], 'k-.')
plt.loglog(runtimes_bench[22], jparerr_bench[22], 'k-')
plt.loglog(runtimes_bench[5], jparerr_bench[5], 'k--')
plt.xlim([4e-3, 2e0])
plt.ylim([1e-14, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta J_{\parallel} $')

plt.subplot(1,2,2)
plt.loglog(evals_bench[15], jparerr_bench[15], 'k-.')
plt.loglog(evals_bench[22], jparerr_bench[22], 'k-')
plt.loglog(evals_bench[5], jparerr_bench[5], 'k--')
plt.xlim([1e4, 1e7])
plt.ylim([1e-14, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta J_{\parallel} $')
plt.tight_layout()

plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[15], Herr_bench[15], 'k-.')
plt.loglog(runtimes_bench[22], Herr_bench[22], 'k-')
plt.loglog(runtimes_bench[5], Herr_bench[5], 'k--')
plt.xlim([4e-3, 2e0])
plt.ylim([1e-14, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta H$')

plt.subplot(1,2,2)
plt.loglog(evals_bench[15], Herr_bench[15], 'k-.')
plt.loglog(evals_bench[22], Herr_bench[22], 'k-')
plt.loglog(evals_bench[5], Herr_bench[5], 'k--')
plt.xlim([1e4, 1e7])
plt.ylim([1e-14, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')
plt.tight_layout()
#%%
# Stellarator
#plt.figure(figsize=(7,3))
#plt.subplot(1,2,1)
#plt.semilogx([0.5, 10], 7.24e-04*np.ones(2), 'tab:gray')
#plt.semilogx(runtimes_bench[1], jparerr_bench[1], '-')
#plt.semilogx(runtimes_bench[3], jparerr_bench[3], ':')
#plt.semilogx(runtimes_bench[4], jparerr_bench[4], '--')
#plt.semilogx(runtimes_bench[5], jparerr_bench[5], '-.')
#plt.semilogx(runtimes_bench[21], jparerr_bench[21], '-.')
#plt.ylim([7.1e-4, 8e-4])
#plt.xlabel('runtime / s')
#plt.ylabel('$\delta J_{\parallel} $')
#plt.subplot(1,2,2)
#plt.semilogx([3e5, 5e6], 7.24e-04*np.ones(2), 'tab:gray')
#plt.semilogx(evals_bench[1], jparerr_bench[1], '-')
#plt.semilogx(evals_bench[3], jparerr_bench[3], ':')
#plt.semilogx(evals_bench[4], jparerr_bench[4], '--')
#plt.semilogx(evals_bench[5], jparerr_bench[5], '-.')
#plt.semilogx(evals_bench[21], jparerr_bench[21], '-.')
#plt.ylim([7.1e-4, 8e-4])
#plt.xlabel('evaluations')
#plt.ylabel('$\delta J_{\parallel} $')
#plt.tight_layout()

#plt.figure(figsize=(7,3))
#plt.subplot(1,2,1)
#plt.loglog(runtimes_bench[1], Herr_bench[1], '-')
#plt.loglog(runtimes_bench[3], Herr_bench[3], ':')
#plt.loglog(runtimes_bench[4], Herr_bench[4], '--')
#plt.loglog(runtimes_bench[5], Herr_bench[5], '-.')
#plt.loglog(runtimes_bench[21], Herr_bench[21], '-.')
#plt.xlabel('runtime / s')
#plt.ylabel('$\delta H$')
#plt.subplot(1,2,2)
#plt.loglog(evals_bench[1], Herr_bench[1], '-')
#plt.loglog(evals_bench[3], Herr_bench[3], ':')
#plt.loglog(evals_bench[4], Herr_bench[4], '--')
#plt.loglog(evals_bench[5], Herr_bench[5], '-.')
#plt.loglog(evals_bench[21], Herr_bench[21], '-.')
#plt.xlabel('evaluations')
#plt.ylabel('$\delta H$')
#plt.tight_layout()
    
#bench.cleanup_bench()
#plt.show()
