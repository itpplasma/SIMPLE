"""
Created: Fri Aug  9 15:50:40 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
#%%
#cmd = '/home/calbert/build/NEO-ORB/bench.x'
#cmd = '/Users/ert/build/NEO-ORB/bench.x'
cmd = '/home/ert/build/NEO-ORB/bench.x'
#cmd = '/bin/echo'

runtimes = []
evalss = []
Herr = []
jparerr = []
runtime = None
evals = None
#%%

def do_run(integ_mode, npoiper2):
    global data, runtimes, evalss, Herr, jparerr, jparmean, Hmean, runtime, evals
    p = Popen(cmd, stdin=PIPE, stdout=PIPE, shell=True)
    
    intxt = r'''&bench
      multi = .False.,
      quasi = .False.,
      tok = .True.,
      integ_mode = {},
      npoiper2 = {},
      nlag = 0,
      ncut = 1000,
      infile = '../RUN/wout_23_1900_fix_bdry.nc',
      outfile = '/tmp/out.txt'
    /          
    '''.format(integ_mode, npoiper2).encode('utf-8')
        
    output = p.communicate(intxt)
    
    rundata = output[0].split(b'\n')[-2].split()
    runtime = float(rundata[0])
    evals = int(rundata[1])
    
    data = np.loadtxt('/tmp/out.txt')
    Hmean = np.mean(data[1:,4])
    jparmean = np.mean(data[1:,5])
    
    runtimes.append(runtime)
    evalss.append(evals)
    Herr.append(np.sqrt(np.mean((data[1:,4]/Hmean-1.0)**2)))
    jparerr.append(np.sqrt(np.mean((data[1:,5]/jparmean-1.0)**2)))

def do_plot():
    plt.figure()
    plt.plot(data[:,0]*np.cos(data[:,1]), data[:,0]*np.sin(data[:,1]), ',')
    
    plt.figure()
    plt.plot(data[1:,4]/Hmean, '-')
    
    plt.figure()
    plt.plot(data[1:,5]/jparmean, '-')
    plt.show()
    
plt.figure()
    
integ_mode = 1
for npoiper2 in [10, 16, 32, 64, 128]:
    do_run(integ_mode, npoiper2)
    #do_plot()
    print('{:.2e} {} {:.2e} {:.2e}'.format(runtime, evals, Herr[-1], jparerr[-1]))
    
plt.loglog(runtimes, Herr, 'x')

runtimes = []
Herr = []
jparerr = []
    
integ_mode = 4
for npoiper2 in [10, 16, 32, 64, 128]:
    do_run(integ_mode, npoiper2)
    #do_plot()
    print('{:.2e} {} {:.2e} {:.2e}'.format(runtime, evals, Herr[-1], jparerr[-1]))
    
plt.loglog(runtimes, Herr, 'o')