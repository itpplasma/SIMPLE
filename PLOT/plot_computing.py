#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 09:56:58 2019

@author: calbert
"""

import numpy as np
import matplotlib.pyplot as plt

cores = np.array([1, 2, 4, 8, 16, 32, 64])
times = np.array([12371.45, 22974.33, 22654.08, 23083.88, 22474.06,
                  22319.14, 25479.76])

speedup = cores/times*1.0/(cores[0]/times[0])

figure()
plt.loglog(cores, speedup, 'k-s')
plt.loglog([1,100],[1,100],'k--')
plt.loglog([1,100],[speedup[1]/2, 50*speedup[1]],'k-.')
plt.xlabel('threads')
plt.ylabel('speedup')
plt.xlim([9e-1,1e2])
plt.ylim([9e-1,1e2])
plt.legend(['actual speedup', 'linear 1', 'linear 2'])

#%% Configuration: multharm=3, tracing time: 2e-3 seconds.

steps_sympl = 2**np.arange(5,12+1)
calls_sympl = np.array([7773841, 13369974, 24780079, 41913639, 69070117, 
                        137083774, 271256655, 541541239])

#calls_rk = np.array([6619030, 9973427, 19374320, 30897192,
#                    45487460, 65510699, 94866328, 140939805,
#                    210298359, 317872723, 486401412])

calls_rk = np.array([30897192, 65510699, 140939805, 317872723])


plt.figure()
plt.semilogy(calls_sympl)
plt.semilogy(calls_rk)

