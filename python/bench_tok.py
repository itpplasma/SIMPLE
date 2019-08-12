"""
Created: Fri Aug  9 15:50:40 2019
@author: Christopher Albert <albert@alumni.tugraz.at>
"""

import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
#%%
cmd = '/home/calbert/build/NEO-ORB/bench.x'

#cmd = '/bin/echo'

p = Popen(cmd, stdin=PIPE, stdout=PIPE, shell=True)
output = p.communicate(b'''&bench
  multi = .False.,
  quasi = .False.,
  integ_mode = 5,
  npoiper2 = 256,
  nlag = 0,
  ncut = 1000,
  outname = '/tmp/out.txt'
/          
''')

#print(output)
rundata = output[0].split(b'\n')[-2].split()
runtime = float(rundata[0])
evals = int(rundata[1])

print(runtime, evals)

#%% Plot results

data = np.loadtxt('/tmp/out.txt')
plt.figure()
plt.plot(data[:,0]*np.cos(data[:,1]), data[:,0]*np.sin(data[:,1]), ',')

plt.figure()
jparmean = np.mean(data[1:,5])
plt.plot(data[1:,5]/jparmean, '-')