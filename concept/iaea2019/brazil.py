from numpy import *
from matplotlib.pyplot import *
from os import path

import exportfig

basedir = path.join(
    '/home/calbert/net/itpplasma/ALPHAPART/SIMPLE-0.9/CLASSPLOTS/',
    'QHs06_10000_new'
)


tpbou = loadtxt(path.join(basedir, 'fort.20000'))
notrace = loadtxt(path.join(basedir, 'fort.10000'))
rapid_p = loadtxt(path.join(basedir, 'fort.10001'))
rapid_t = loadtxt(path.join(basedir, 'fort.10002'))
regul_p = loadtxt(path.join(basedir, 'fort.10011'))
regul_t = loadtxt(path.join(basedir, 'fort.10012'))
stoch_p = loadtxt(path.join(basedir, 'fort.10021'))
stoch_t = loadtxt(path.join(basedir, 'fort.10022'))

nskip = 3
nskip2 = 1
phase = 0
regul_p[:, 0] = mod(regul_p[:, 0] + phase, 2*pi) - phase
regul_t[:, 0] = mod(regul_t[:, 0] + phase, 2*pi) - phase
rapid_p[:, 0] = mod(rapid_p[:, 0] + phase, 2*pi) - phase
rapid_t[:, 0] = mod(rapid_t[:, 0] + phase, 2*pi) - phase
tpbou[:, 0] = mod(tpbou[:, 0] + phase, 2*pi) - phase
#stoch_p[:, 0] = mod(stoch_p[:, 0] + phase, 2*pi)
stoch_t[:, 0] = mod(stoch_t[:, 0] + phase, 2*pi) - phase

figure(figsize=(1.9, 1.7))
scatter(regul_p[::nskip, 0]/pi, regul_p[::nskip, 1], marker='s', color='#009C3B')
scatter(regul_t[::nskip, 0]/pi, regul_t[::nskip, 1], marker='s', color='#009C3B')
scatter(rapid_p[::nskip, 0]/pi, rapid_p[::nskip, 1], s=30, marker='o',
        facecolors='none', edgecolors='#FFDF00')
scatter(rapid_t[::nskip, 0]/pi, rapid_t[::nskip, 1], s=30, marker='o',
        facecolors='none', edgecolors='#FFDF00')
plot(tpbou[::nskip, 0]/pi, tpbou[::nskip, 1], 'ws', markersize=0.5)
plot(tpbou[::nskip, 0]/pi, -tpbou[::nskip, 1], 'ws', markersize=0.5)
#scatter(stoch_p[:, 0]/pi, stoch_p[:, 1], marker='x', color='#002776')
scatter(stoch_t[::nskip2, 0]/pi, stoch_t[::nskip2, 1], marker='x', color='#002776')
#xlim([-1,1])
xlim([0,2])
ylim([-0.8, 0.8])

xlabel(r'$\vartheta / \pi$', labelpad=0)
ylabel(r'$v_\parallel / v$', labelpad=-4)

exportfig.exporteps('brazil_qh_s06')
