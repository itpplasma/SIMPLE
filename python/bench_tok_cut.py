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

bench.multi = 0
bench.quasi = 0
bench.tok = 1
bench.integ_mode = 1
bench.npoiper2 = 16
bench.nlag = 0
bench.nplagr_invar = 4
bench.ncut = 10000
bench.nt = 10000
bench.rtol = 1e-13
bench.init_bench()
    
#%% Benchmarking parallel invariant
bench.ncut = 10000

# RK45
bench.quasi = 1
bench_multi = 0
tols = np.array([1e-6, 1e-8, 1e-10, 1e-12, 1e-14])
npoipers0 = np.array([5, 8, 16, 32, 64])
#npoipers = np.array([8, 16, 32, 64, 128])
nplagr_invar = [2, 2, 4, 6, 8]
integ_mode = 0
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers0, nplagr_invar, tols)
#%%
bench.quasi = 0
bench.rtol = 1e-13

bench.multi = 1
# Stellarator
#npoipers = [20, 24, 28, 32, 48, 64, 128, 256, 512]

# Tokamak
bench.nt = 10000

#npoipers = [32, 64, 128, 256, 512]
#nplagr_invar = np.concatenate([[6], 6*np.ones(len(npoipers)-1, int)])
#integ_mode = 22
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

npoipers = [8, 16, 32, 64, 128, 256]
nplagr_invar = np.concatenate([[2,4,6], 6*np.ones(len(npoipers)-3, int)])
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

npoipers = [5, 8, 16, 32, 64, 128, 256]
nplagr_invar = np.concatenate([[2,2,4,6], 6*np.ones(len(npoipers)-3, int)])
integ_mode = 3
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 4
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)
#
#integ_mode = 15
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

# Stellarator
#npoipers = [20, 24, 28, 32, 64, 128]
# Tokamak
#npoipers = [5, 8, 16, 32, 64, 128, 256]
#nplagr_invar = np.concatenate([[2,2,4,6], 6*np.ones(len(npoipers)-3, int)])
#integ_mode = 5
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

#integ_mode = 6
#(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
#  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)
    

#%%
# Tokamak
from exportfig import *

plt.figure()
plt.loglog(runtimes_bench[0]/bench.ncut, jparerr_bench[0], '-', color='tab:gray')
plt.loglog(runtimes_bench[1]/bench.ncut, jparerr_bench[1], 'k-')
plt.loglog(runtimes_bench[3]/bench.ncut, jparerr_bench[3], 'k:')
plt.loglog(runtimes_bench[4]/bench.ncut, jparerr_bench[4], 'k--')
plt.loglog(runtimes_bench[21]/bench.ncut, jparerr_bench[21], 'k-.')
plt.xlim([3e-6, 1e-3])
plt.ylim([1e-13, 1e-1])
plt.xlabel('CPU time per bounce period / s')
plt.ylabel(r'$\delta J_{\parallel}  / \bar{J_{\parallel}}$')

exportpng('fig_tok_cut1')


plt.figure()
plt.loglog(evals_bench[0]/bench.ncut, jparerr_bench[0], '-', color='tab:gray')
plt.loglog(evals_bench[1]/bench.ncut, jparerr_bench[1], 'k-')
plt.loglog(evals_bench[3]/bench.ncut, jparerr_bench[3], 'k:')
plt.loglog(evals_bench[4]/bench.ncut, jparerr_bench[4], 'k--')
plt.loglog(evals_bench[21]/bench.ncut, jparerr_bench[21], 'k-.')
#plt.xlim([1e4, 1e7])
plt.xlim([1e1, 1e4])
plt.ylim([1e-13, 1e-1])
plt.xlabel('evaluations $N_{\mathrm{f}}$ per bounce period')
plt.ylabel(r'$\delta J_{\parallel}  / \bar{J_{\parallel}}$')

exportpng('fig_tok_cut2')

#%%

plt.figure()
plt.loglog(runtimes_bench[0]/bench.ncut, Herr_bench[0], '-', color='tab:gray')
plt.loglog(runtimes_bench[1]/bench.ncut, Herr_bench[1], 'k-')
plt.loglog(runtimes_bench[3]/bench.ncut, Herr_bench[3], 'k:')
plt.loglog(runtimes_bench[4]/bench.ncut, Herr_bench[4], 'k--')
plt.loglog(runtimes_bench[21]/bench.ncut, Herr_bench[21], 'k-.')
#plt.xlim([4e-3, 2e0])
plt.ylim([1e-14, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta H$')

plt.figure()
plt.loglog(evals_bench[0]/bench.ncut, Herr_bench[0], '-', color='tab:gray')
plt.loglog(evals_bench[1]/bench.ncut, Herr_bench[1], 'k-')
plt.loglog(evals_bench[3]/bench.ncut, Herr_bench[3], 'k:')
plt.loglog(evals_bench[4]/bench.ncut, Herr_bench[4], 'k--')
plt.loglog(evals_bench[21]/bench.ncut, Herr_bench[21], 'k-.')
#plt.xlim([1e4, 1e7])
plt.ylim([1e-14, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')
#    
#plt.figure(figsize=(7,3))
#plt.subplot(1,2,1)
#plt.loglog(runtimes_bench[15]/bench.ncut, jparerr_bench[15], 'k-.')
#plt.loglog(runtimes_bench[22]/bench.ncut, jparerr_bench[22], 'k-')
#plt.loglog(runtimes_bench[5]/bench.ncut, jparerr_bench[5], 'k--')
#plt.xlim([4e-3, 2e0])
#plt.ylim([1e-14, 1e-1])
#plt.xlabel('runtime / s')
#plt.ylabel('$\delta J_{\parallel} $')
#
#plt.subplot(1,2,2)
#plt.loglog(evals_bench[15]/bench.ncut, jparerr_bench[15], 'k-.')
#plt.loglog(evals_bench[22]/bench.ncut, jparerr_bench[22], 'k-')
#plt.loglog(evals_bench[5]/bench.ncut, jparerr_bench[5], 'k--')
#plt.xlim([1e4, 1e7])
#plt.ylim([1e-14, 1e-1])
#plt.xlabel('evaluations')
#plt.ylabel('$\delta J_{\parallel} $')
#plt.tight_layout()
#
#plt.figure(figsize=(7,3))
#plt.subplot(1,2,1)
#plt.loglog(runtimes_bench[15]/bench.ncut, Herr_bench[15], 'k-.')
#plt.loglog(runtimes_bench[22]/bench.ncut, Herr_bench[22], 'k-')
#plt.loglog(runtimes_bench[5]/bench.ncut, Herr_bench[5], 'k--')
#plt.xlim([4e-3, 2e0])
#plt.ylim([1e-14, 1e-1])
#plt.xlabel('runtime / s')
#plt.ylabel('$\delta H$')
#
#plt.subplot(1,2,2)
#plt.loglog(evals_bench[15]/bench.ncut, Herr_bench[15], 'k-.')
#plt.loglog(evals_bench[22]/bench.ncut, Herr_bench[22], 'k-')
#plt.loglog(evals_bench[5]/bench.ncut, Herr_bench[5], 'k--')
#plt.xlim([1e4, 1e7])
#plt.ylim([1e-14, 1e-1])
#plt.xlabel('evaluations')
#plt.ylabel('$\delta H$')
#plt.tight_layout()


#%%
