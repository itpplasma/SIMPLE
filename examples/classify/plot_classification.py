import numpy as np
import matplotlib.pyplot as plt
import os

def doplot_inner(prompt, regular, stochastic, bminmax, relvar_jpar, outfile,
    arrows):
    tot1=prompt+regular+stochastic
    tot2=regular+stochastic
    wpr=-1.
    wst=0.
    wre=1.
    cla1=wre*regular+wst*stochastic+wpr*prompt
    cla1[tot1==0]=np.nan
    tot1[tot1==0]=np.nan
    cla1=cla1/tot1

    relvar_jpar[tot2==0]=np.nan
    tot2[tot2==0]=np.nan
    jpvar=relvar_jpar/tot2

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
    for a in arrows:
        plt.arrow(a[0], a[1], a[2], a[3], length_includes_head=True,
            head_width=0.01, head_length=0.04, ec='w', fc='w')
    plt.xticks([0,.25,.5,.75,1])
    plt.ylim(jpmin, 1)
    plt.xlabel('$s$')
    plt.ylabel('$J_\perp$')

    plt.tight_layout()
    plt.savefig(outfile)

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
    relvar_jpar = np.loadtxt(os.path.join(basedir, 'relvar_jpar.dat'))

    doplot_inner(prompt1, regular1, stochastic1, bminmax, relvar_jpar, outfile1, [])
    doplot_inner(prompt2, regular2, stochastic2, bminmax, relvar_jpar, outfile2,arrows)


doplot('../fig/classification_new/QI_Subbotin',
    '../../paper/fig/qi_S_full_jp.pdf',
    '../../paper/fig/qi_S_full_id.pdf',
    [(.25, 0.875, 0.45-0.25, 0.0), (.25, 0.56, 1.0-0.25, 0.0)])

doplot('../fig/classification_new/QI_Drevlak',
    '../../paper/fig/qi_D_full_jp.pdf',
    '../../paper/fig/qi_D_full_id.pdf', [])

doplot('../fig/classification_new/QH_Drevlak/',
    '../../paper/fig/qh_D_full_jp.pdf',
    '../../paper/fig/qh_D_full_id.pdf', [])

doplot('../fig/classification_new/QA_Henneberg/3.5Mev',
    '../../paper/fig/qa_H_full_jp.pdf',
    '../../paper/fig/qa_H_full_id.pdf', [])

doplot('../fig/classification_new/QA_Henneberg/35Kev',
    '../../paper/fig/qa_H_red_jp.pdf',
    '../../paper/fig/qa_H_red_id.pdf', [])

doplot('../fig/classification_new/QA_Landreman/LOWRES',
    '../../paper/fig/qa_L_full_jp.pdf',
    '../../paper/fig/qa_L_full_id.pdf', [])
