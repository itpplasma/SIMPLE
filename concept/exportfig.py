import matplotlib.pyplot as plt
import os

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[3.2, 2.4])
plt.rc('lines', linewidth=0.5)
plt.rc('lines', markersize=3)

def exportfig(filename):
    print('export: ' + filename)
    plt.tight_layout()
    plt.savefig(filename + '.pdf')
    os.system('pdfcrop ' + filename + '.pdf ' + filename + '.pdf')
    

def exportpng(filename):
    print('export: ' + filename)
    plt.tight_layout()
    plt.savefig(filename + '.png', dpi=300)
    os.system(r'mogrify -trim ' + filename + '.png ')
    