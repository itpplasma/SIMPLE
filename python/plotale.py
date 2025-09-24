from numpy import cos, sin, empty_like, linspace, roll
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure, plot, xlabel, ylabel, subplot, grid, legend



def plot_orbit(z, ax):

    ax.plot(z[0, :] * cos(z[1, :]),
            z[0, :] * sin(z[1, :]),
            'o', markersize=1)

    ax.set_xlabel(r'$R-R_0$')
    ax.set_ylabel(r'$Z$')
    return ax


def plot_mani(z):
    # Torus parameters
    R0 = 1  # Distance from the center of the tube to the center of the torus
    r = 0.5  # Radius of the tube

    # Create meshgrid
    theta = linspace(0, 2 * np.pi, 100)
    phi = linspace(0.3, 1.4*np.pi, 100)
    theta, phi = np.meshgrid(theta, phi)

    # Parametric equations for a torus
    X = (R0 + r * cos(theta)) * cos(phi)
    Y = (R0 + r * cos(theta)) * sin(phi)
    Z = r * sin(theta)

    # Plotting
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(X, Y, Z, cmap='plasma', edgecolor='none')

    x = (R0 + z[0,:]*cos(z[1,:]))*cos(z[2,:])
    y = (R0 + z[0,:]*cos(z[1,:]))*sin(z[2,:])
    z = z[0,:]*sin(z[1,:])
    ax.scatter(x,y,z)
    ax.set_zlim([-1,1])
    






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
