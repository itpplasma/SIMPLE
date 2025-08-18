"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""
from numpy import cos, sin, empty_like, linspace, roll
from matplotlib.pyplot import figure, plot, xlabel, ylabel, subplot, grid, legend

def plot_orbit(z):
    figure()
    plot(z[0,:]*cos(z[1,:]), z[0,:]*sin(z[1,:]),',')
    xlabel(r'$R-R_0$')
    ylabel(r'$Z$')


def plot_cost_function(F, zlast, zold, pthold):
    """Plot the cost function"""
    r = linspace(-0.2,0.2,100)
    fres = empty_like(r)
    for kr in range(len(r)):
        fres[kr] = F([r[kr]], zold[1:], pthold)
    
    figure()
    subplot(2,1,1)
    plot(r,fres)
    out_conv = F([zlast[0]], zold[1:], pthold)
    plot(zlast[0], out_conv, 'rx')
    xlabel('r')
    ylabel('F')
    grid()

def plot_cost_function_jac(F, zlast, zold, pthold):
    """Plot the cost function with Jacobian"""
    
    r = linspace(-0.2,0.2,100)
    fres = empty_like(r)
    jac = empty_like(r)
    for kr in range(len(r)):
        out = F([r[kr]], zold[1:], pthold)
        fres[kr] = out[0]
        jac[kr] = out[1][0]
        
    jacnum = (roll(fres,-1)-roll(fres,1))/(2.*(r[1]-r[0]))
    jacnum[0] = jacnum[1]
    jacnum[-1] = jacnum[-2]
    
    figure()
    subplot(2,1,1)
    plot(r,fres)
    out_conv = F([zlast[0]],zold[1:],pthold)
    plot(zlast[0],out_conv[0],'rx')
    xlabel('r')
    ylabel('F')
    grid()
    
    subplot(2,1,2)
    plot(r,jac)
    plot(r,jacnum,'r--')
    plot(zlast[0],out_conv[1],'rx')
    xlabel('r')
    ylabel('dF/dr')
    legend(['analytical', 'numerical'])
    grid()
