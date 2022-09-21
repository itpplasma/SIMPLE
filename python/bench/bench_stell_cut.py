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
bench.tok = 0
bench.integ_mode = 1
bench.npoiper2 = 16
bench.nlag = 0
bench.nplagr_invar = 4
bench.ncut = 10000
bench.nt = 10000
bench.rtol = 1e-13
bench.init_bench()
#%%
bench.ncut = 10000

# RK45
bench.quasi = 1
bench.npoiper = 16
bench_multi = 0
tols = np.array([1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 1e-10])
npoipers0 = np.array([16, 16, 16, 32, 32, 64])
nplagr_invar = [4, 4, 4, 6, 6, 8]
integ_mode = 0
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers0, nplagr_invar, tols)
#%%
bench.quasi = 0
bench.rtol = 1e-13

bench.multi = 1
# Stellarator
npoipers = [19, 20, 21, 22, 23, 24, 26, 28, 32, 48, 64]#, 128, 256, 512]

nplagr_invar = np.concatenate([4*np.ones(len(npoipers), int)])

integ_mode = 21
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

bench.multi = 0

integ_mode = 1
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 3
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

integ_mode = 4
(runtimes_bench[integ_mode], evals_bench[integ_mode], Herr_bench[integ_mode],
  jparerr_bench[integ_mode]) = do_cutbench(integ_mode, npoipers, nplagr_invar)

#%%
# Stellarator
from exportfig import *

plt.figure()
plt.semilogx([1e-4, 3e-3], 7.24e-04*np.ones(2), 'lightgray')
plt.loglog(runtimes_bench[0]/bench.ncut, jparerr_bench[0], '-', color='gray')
plt.loglog(runtimes_bench[1][1:]/bench.ncut, jparerr_bench[1][1:], 'k-')
plt.loglog(runtimes_bench[3][1:]/bench.ncut, jparerr_bench[3][1:], 'k:')
plt.loglog(runtimes_bench[4]/bench.ncut, jparerr_bench[4], 'k--')
#plt.semilogx(runtimes_bench[5], jparerr_bench[5], '-.')
plt.loglog(runtimes_bench[21][1:]/bench.ncut, jparerr_bench[21][1:], 'k-.')
plt.xlim([1e-4, 3e-3])
plt.ylim([6e-4, 1e-1])
plt.xlabel('CPU time per bounce period / s')
plt.ylabel(r'$\delta J_{\parallel}  / \bar{J_{\parallel}}$')
exportpng('fig_stell_cut1')

plt.figure()
plt.loglog([8e1, 2e3], 7.24e-04*np.ones(2), 'lightgray')
plt.loglog(evals_bench[0]/bench.ncut, jparerr_bench[0], '-', color='gray')
plt.loglog(evals_bench[1][1:]/bench.ncut, jparerr_bench[1][1:], 'k-')
plt.loglog(evals_bench[3][1:]/bench.ncut, jparerr_bench[3][1:], 'k:')
plt.loglog(evals_bench[4]/bench.ncut, jparerr_bench[4], 'k--')
#plt.semilogx(evals_bench[5], jparerr_bench[5], '-.')
plt.loglog(evals_bench[21][1:]/bench.ncut, jparerr_bench[21][1:], 'k-.')
plt.xlim([8e1, 2e3])
plt.ylim([6e-4, 1e-1])
plt.xlabel('evaluations $N_{\mathrm{f}}$ per bounce period')
plt.ylabel(r'$\delta J_{\parallel}  / \bar{J_{\parallel}}$')
exportpng('fig_stell_cut2')

#%%
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[1]/bench.ncut, Herr_bench[1], '-')
plt.loglog(runtimes_bench[3]/bench.ncut, Herr_bench[3], ':')
plt.loglog(runtimes_bench[4]/bench.ncut, Herr_bench[4], '--')
#plt.loglog(runtimes_bench[5], Herr_bench[5], '-.')
plt.loglog(runtimes_bench[21]/bench.ncut, Herr_bench[21], '-.')
plt.xlabel('runtime / s')
plt.ylabel('$\delta H$')
plt.subplot(1,2,2)
plt.loglog(evals_bench[1]/bench.ncut, Herr_bench[1], '-')
plt.loglog(evals_bench[3]/bench.ncut, Herr_bench[3], ':')
plt.loglog(evals_bench[4]/bench.ncut, Herr_bench[4], '--')
#plt.loglog(evals_bench[5], Herr_bench[5], '-.')
plt.loglog(evals_bench[21]/bench.ncut, Herr_bench[21], '-.')
plt.xlabel('evaluations')
plt.ylabel('$\delta H$')
plt.tight_layout()
    


#%%
