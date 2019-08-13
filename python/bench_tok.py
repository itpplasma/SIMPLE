"""
Created: Fri Aug  9 15:50:40 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from fffi import fortran_library, fortran_module

runtimes = []
evalss = []
Herr = []
jparerr = []
runtimes_bench = {}
evals_bench = {}
Herr_bench = {}
jparerr_bench = {}

runtime = None
evals = None

libneo_orb = fortran_library('neo_orb', path='../lib')

bench = fortran_module(libneo_orb, 'neo_orb_bench')
bench.fdef("""  
    integer :: npoiper2
    double precision :: rbig, dtau, dtaumax
    
    integer :: nt
    double precision, allocatable :: out(:, :)
    
    double precision :: starttime, endtime
    
    
    logical :: multi
    logical :: quasi
    logical :: tok
    integer :: integ_mode
    integer :: nlag
    integer :: nplagr_invar
    integer :: ncut
    
    double precision, allocatable :: var_cut(:, :)
    
    double precision :: taub
    
    
    subroutine init_bench
    end
    
    subroutine do_bench
    end
    
    subroutine cleanup_bench
    end
    
    subroutine test_cuts(nplagr)
      integer, intent(in) :: nplagr
    end
""")

diag = fortran_module(libneo_orb, 'diag_mod')
diag.fdef("""  
    integer :: icounter
""")

libneo_orb.compile()
bench.load()
diag.load()
#%%
bench.multi = 0
bench.quasi = 0
bench.tok = 1
bench.integ_mode = 1
bench.npoiper2 = 64
bench.nlag = 0
bench.nplagr_invar = 4
bench.ncut = 1000
bench.init_bench()

#%%

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
        if(krun == 1): do_plot()
        print('{} {:.2e} {} {:.2e} {:.2e}'.format(npoipers[krun],
            runtime, evals, Herr[-1], jparerr[-1]))
    
    return runtimes, evalss, Herr, jparerr
    
#%%
bench.multi = 1
# Stellarator
#npoipers = [20, 24, 28, 32, 48, 64, 128, 256, 512]

# Tokamak

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


#bench.cleanup_bench()
    
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
    
# Tokamak
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.loglog(runtimes_bench[1], jparerr_bench[1], 'k-')
plt.loglog(runtimes_bench[3], jparerr_bench[3], 'k:')
plt.loglog(runtimes_bench[4], jparerr_bench[4], 'k--')
plt.loglog(runtimes_bench[21], jparerr_bench[21], 'k-.')
plt.xlim([4e-3, 2e0])
plt.ylim([1e-14, 1e-1])
plt.xlabel('runtime / s')
plt.ylabel('$\delta J_{\parallel} $')

plt.subplot(1,2,2)
plt.loglog(evals_bench[1], jparerr_bench[1], 'k-')
plt.loglog(evals_bench[3], jparerr_bench[3], 'k:')
plt.loglog(evals_bench[4], jparerr_bench[4], 'k--')
plt.loglog(evals_bench[21], jparerr_bench[21], 'k-.')
plt.xlim([1e4, 1e7])
plt.ylim([1e-14, 1e-1])
plt.xlabel('evaluations')
plt.ylabel('$\delta J_{\parallel} $')
plt.tight_layout()

plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
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
    