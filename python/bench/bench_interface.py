"""
Created: Tue Aug 20 15:07:28 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
from fffi import fortran_library, fortran_module

libneo_orb = fortran_library('neo_orb', path='../lib')

bench = fortran_module(libneo_orb, 'simple_bench')
bench.fdef("""
    integer :: npoiper2
    double precision :: rbig, dtau, dtaumax

    integer :: nt

    double precision :: starttime, endtime


    logical :: multi
    logical :: quasi
    logical :: tok
    integer :: integ_mode
    integer :: nlag
    integer :: nplagr_invar
    integer :: ncut

    double precision :: taub
    double precision :: rtol


    subroutine init_bench
    end

    subroutine do_bench
    end

    subroutine cleanup_bench
    end

    subroutine test_cuts(nplagr)
      integer, intent(in) :: nplagr
    end

    subroutine minsqdist(za, zref, result)
      double precision, intent(in) :: za(:,:)
      double precision, intent(in) :: zref(:,:)
      double precision, intent(out) :: result(:)
    end
""")

diag = fortran_module(libneo_orb, 'diag_mod')
diag.fdef("""
    integer :: icounter
""")

libneo_orb.compile()
bench.load()
diag.load()

def do_run(integ_mode, npoiper2, nplagr_invar, rtol = None):

    bench.integ_mode = integ_mode
    bench.npoiper2 = npoiper2
    if rtol is not None: bench.rtol = rtol
    bench.nplagr_invar = nplagr_invar

    bench.do_bench()

    runtime = bench.endtime - bench.starttime
    evals = diag.icounter

    data = np.loadtxt('/tmp/out.txt')
    Hmean = np.sqrt(np.mean(data[1:,4]**2))
    jparmean = np.sqrt(np.mean(data[1:,5]**2))

    Herr = np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2))
    jparerr = np.sqrt(np.mean((data[1:,5]/jparmean-1.0)**2))

    return runtime, evals, Herr, jparerr


def do_cutbench(integ_mode, npoipers, nplagr_invar, rtols = None):
    runtimes = []
    evalss = []
    Herrs = []
    jparerrs = []
    nrun = len(npoipers)
    print('')
    for krun in range(nrun):
        if rtols is not None:
            rtol = rtols[krun]
        else:
            rtol = None

        runtime, evals, Herr, jparerr = do_run(
                integ_mode, npoipers[krun], nplagr_invar[krun], rtol)
        runtimes.append(runtime)
        evalss.append(evals)
        Herrs.append(Herr)
        jparerrs.append(jparerr)
        print('{} {:.2e} {} {:.2e} {:.2e}'.format(npoipers[krun],
            runtime, evals, Herr, jparerr))

    return np.array(runtimes), np.array(evalss), np.array(Herrs), np.array(jparerrs)


#%%

def do_singlebench(integ_mode, npoipers, nskip, zref):
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
        avgdists.append(avgdist(np.stack((x,y)).copy(order='F'), zref))
        Hmean = np.sqrt(np.mean(data[1:,4]**2))
        Herr.append(np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2)))
        runtimes.append(bench.endtime - bench.starttime)
        evalss.append(diag.icounter)

    return np.array(runtimes), np.array(evalss), np.array(avgdists), np.array(Herr)

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

def avgdist(z, zref):
    nz = z.shape[1]
    sqdists = np.empty(nz)
    bench.minsqdist(z, zref, sqdists)
    #print(sqdists)
    return np.sqrt(np.mean(sqdists))
