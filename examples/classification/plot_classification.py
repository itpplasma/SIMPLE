import numpy as np
import matplotlib.pyplot as plt
import os

def doplot_inner(prompt, regular, stochastic, bminmax, outfile, arrows):
    tot1=prompt+regular+stochastic
    wpr=-1.
    wst=0.
    wre=1.
    cla1=wre*regular+wst*stochastic+wpr*prompt
    cla1[tot1==0]=np.nan
    tot1[tot1==0]=np.nan
    cla1=cla1/tot1

    bmin=min(bminmax[:,1])
    bmax=max(bminmax[:,2])
    jpmin=bmin/bmax
    print("Minimum J_perp: ", jpmin)

    plt.figure(figsize=(5,4))
    plt.imshow(cla1.T, origin='lower', extent=[0,1,0,1], aspect='auto',
        vmin=-1.0, vmax=1.0, clim=(-1.0,1.0))
    plt.plot(bminmax[:,0],bmin/bminmax[:,1],'k')
    plt.plot(bminmax[:,0],bmin/bminmax[:,2],'k--')
    plt.plot([.25, .25],[0.0, bmin/bminmax[250,1]], 'k', lw=0.75)
    for a in arrows:
        plt.arrow(a[0], a[1], a[2], a[3], length_includes_head=True,
            head_width=0.01, head_length=0.04, ec='w', fc='w')
    plt.xticks([0,.25,.5,.75,1])
    plt.ylim(jpmin, 1)
    plt.xlabel('Normalized toroidal flux $s$')
    plt.ylabel('Perpendicular invariant $J_\perp$')

    plt.tight_layout()

def doplot(basedir, outfile1, outfile2, arrows):
    # J_parallel regular/chaotic classification
    prompt1 = np.loadtxt(os.path.join(basedir, 'prompt1.dat'))
    regular1 = np.loadtxt(os.path.join(basedir, 'regular1.dat'))
    stochastic1 = np.loadtxt(os.path.join(basedir, 'stochastic1.dat'))

    # Topological ideal/non-ideal classification
    prompt2 = np.loadtxt(os.path.join(basedir, 'prompt2.dat'))
    regular2 = np.loadtxt(os.path.join(basedir, 'regular2.dat'))
    stochastic2 = np.loadtxt(os.path.join(basedir, 'stochastic2.dat'))

    bminmax = np.loadtxt(os.path.join(basedir, 'bminmax.dat'))

    doplot_inner(prompt1, regular1, stochastic1, bminmax, outfile1, [])
    plt.title(r'$J_\parallel$ classifier')
    cb = plt.colorbar()
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(['Prompt losses', 'Non-ideal', 'Prompt'])
    plt.savefig(outfile1)
    doplot_inner(prompt2, regular2, stochastic2, bminmax, outfile2, arrows)
    plt.title('Topological classifier')
    # Set colorbar tick labels
    cb = plt.colorbar()
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(['Prompt losses', 'Non-ideal', 'Ideal'])
    plt.savefig(outfile2)

doplot('.', 'class_jpar.pdf', 'class_ideal.pdf', [])
