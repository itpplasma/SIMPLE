import numpy as np
import matplotlib.pyplot as plt
import os

def doplot(prompt, regular, stochastic, bminmax):
    tot1=prompt+regular+stochastic

    wpr=-1.
    wst=0.
    wre=1.
    cla1=wre*regular+wst*stochastic+wpr*prompt
    cla1[tot1==0]=np.nan
    tot1[tot1==0]=np.nan
    cla1=cla1/tot1

    ns=np.shape(cla1)[0]
    s=(np.arange(ns)+0.5)/ns
    nj=np.shape(cla1)[1]
    bmin=min(bminmax[:,1])
    bmax=max(bminmax[:,2])
    jp=(np.arange(nj)+0.5)/nj
    jpmin=bmin/bmax
    print(jpmin)

    plt.figure(figsize=(3,3))
    plt.imshow(cla1.T, origin='lower', extent=[0,1,0,1], aspect='auto',
        vmin=-1.0, vmax=1.0, clim=(-1.0,1.0))
    plt.plot(bminmax[:,0],bmin/bminmax[:,1],'k')
    plt.plot(bminmax[:,0],bmin/bminmax[:,2],'k--')
    plt.plot([.25, .25],[0.0, bmin/bminmax[250,1]], 'k', lw=0.75)

    plt.xticks([0,.25,.5,.75,1])
    plt.ylim(jpmin, 1)
    plt.xlabel('$s$')
    plt.ylabel('$J_\perp$')

    plt.tight_layout()

# J_parallel regular/chaotic classification
prompt = np.loadtxt('fort.10002')
ideal = np.loadtxt('fort.40012')
nonideal = np.loadtxt('fort.40022')


bminmax = np.loadtxt('bminmax.dat')
# relvar_jpar = np.loadtxt(os.path.join(basedir, 'relvar_jpar.dat'))

doplot(prompt, ideal, nonideal, bminmax)
