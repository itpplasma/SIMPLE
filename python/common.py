"""
Created:  2018-08-08
Modified: 2019-03-07
Author:   Christopher Albert <albert@alumni.tugraz.at>
"""

from numpy import empty, nan, zeros, ones, sqrt, mod, array
from field_test import field
#%%

taub = 7800 # estimated bounce time

def timesteps(steps_per_bounce, nbounce):
    """ compute step size and number of timesteps """
    return taub/steps_per_bounce, nbounce*steps_per_bounce

f = field()

qe = 1.0; m = 1.0; c = 1.0 # particle charge, mass and speed of light
mu = 1e-5 # magnetic moment

# Initial conditions
r0    = 0.1
th0   = 1.5
ph0   = 0.0
vpar0 = 0.0

# Compute toroidal momentum from initial conditions
f.evaluate(r0, th0, ph0)  
pph0 = m*vpar0*f.hph + qe/c*f.Aph

# Buffer to avoid double evaluation within black-box root finding routines
bufsize = 16
kbuf = 0
zbuf = empty((4,bufsize)); zbuf[:] = nan
neval = 0
buffer = empty(bufsize,dtype=object)
fbuf = empty(bufsize,dtype=object)

dHdvpar = lambda vpar: m*vpar
pth = lambda vpar,hth,Ath: m*vpar*hth + qe*Ath/c
pph = lambda vpar,hph,Aph: m*vpar*hph + qe*Aph/c

def get_val(z):
    """ evaluates H, pth, vpar at z """
    f.evaluate(z[0], z[1], z[2])  
    vpar = 1.0/f.hph*(z[3]-qe/c*f.Aph)
    H = m*vpar**2/2 + mu*f.B + qe*f.Phie
    pth = m*f.hth*vpar + qe/c*f.Ath
    ret = [H, pth, vpar]
    return ret
    
def get_der(z):
    """ evaluates H, pth, vpar and first derivatives at z """
    global neval, kbuf
    
    # re-evaluation at same point
    if(bufsize>0):
        for kb in range(bufsize):
            if all(z==zbuf[:,kb]):
                f.__dict__ = fbuf[kb].__dict__.copy()
                return buffer[kb]
        
    neval = neval+1
    
    [H,pth,vpar] = get_val(z)
    
    dvpardx = -(qe/(m*c)*f.dAph + vpar*f.dhph)/f.hph
    dvpardpph = 1.0/(m*f.hph)
        
    dHdx = m*vpar*dvpardx + mu*f.dB + qe*f.dPhie
    dHdpph = m*vpar/f.hph
    
    dpthdx = m*dvpardx*f.hth + m*vpar*f.dhth + qe/c*f.dAth
    dpthdpph = f.hth/f.hph
    
    ret = [H, pth, vpar, dHdx, dHdpph, dpthdx, dpthdpph, dvpardx, dvpardpph]
    
    # fill buffer
    if(bufsize>0):
        zbuf[:,kbuf] = z[:]
        buffer[kbuf] = ret
        fbuf[kbuf] = field()
        fbuf[kbuf].__dict__ = f.__dict__.copy()
        kbuf = mod(kbuf+1,bufsize)
        
    return ret

def get_der2(z):
    """ evaluates H, pth, vpar with first and second derivatives at z """
    [H, pth, vpar, 
     dHdx, dHdpph, dpthdx, dpthdpph, dvpardx, dvpardpph] = get_der(z)
    
    d2pthdx2 = zeros(6)    # d2dr2, d2dth2, d2dph2, d2drdth, d2drdph, d2dthdph
    d2pthdpphdz = zeros(4) # d2dpphdr, d2dpphdth, d2dpphdph, d2dpph2
    d2Hdx2 = zeros(6)
    d2Hdpphdz = zeros(4)
    d2vpardx2 = zeros(6)
    d2vpardpphdz = zeros(4)
    
    d2vpardx2[:3] = -(qe/(m*c)*f.d2Aph[:3] + f.d2hph[:3]*vpar 
                      + 2*f.dhph*dvpardx)/f.hph
    d2vpardx2[3]  = -(qe/(m*c)*f.d2Aph[3] + f.d2hph[3]*vpar 
                      + f.dhph[0]*dvpardx[1] + f.dhph[1]*dvpardx[0])/f.hph
    d2vpardpphdz[:3] = -1.0/(m*f.hph**2)*f.d2hph[:3]
             
    d2pthdx2[:3] = m*(d2vpardx2[:3]*f.hth + 2*dvpardx*f.dhth + vpar*f.d2hth[:3]
                      + qe/(m*c)*f.d2Ath[:3])
    
    d2pthdx2[3] = m*(d2vpardx2[3]*f.hth + dvpardx[0]*f.dhth[1] 
                     + dvpardx[1]*f.dhth[0] + vpar*f.d2hth[3] + qe/(m*c)*f.d2Ath[3])
    d2pthdpphdz[:3] = f.dhth/f.hph - f.hth/f.hph**2*f.dhph
    
    d2Hdx2[:3] = m*(dvpardx[:3]**2 + vpar*d2vpardx2[:3]) \
                 + mu*f.d2B[:3] + qe*f.d2Phie[:3]
    d2Hdx2[3]  = m*(dvpardx[0]*dvpardx[1] + vpar*d2vpardx2[3]) \
                 + mu*f.d2B[3] + qe*f.d2Phie[3]
    d2Hdpphdz[:3] = m*(1.0/f.hph*dvpardx - vpar/f.hph**2*f.dhph)
    
    return [H, pth, vpar, 
            dHdx, dHdpph, dpthdx, dpthdpph, dvpardx, dvpardpph,
            d2pthdx2, d2pthdpphdz, d2Hdx2, d2Hdpphdz, d2vpardx2, d2vpardpphdz]
    
class NewtonResult():
    x = array([])
    xold = array([])
    
def newton1(f, x0, rtol, atol, args):
    ret = NewtonResult()
    x = array([x0])
    xold = array([x0*(1e0+1e30*rtol)])
    fval = atol*array([1e30,1e30])
    while(abs(fval[0]) > atol and abs(x-xold)/(abs(x)*(1e0+1e-30)) > rtol):
        fval = f(x, *args)
        xold = x
        x = x - fval[0]/fval[1]
    ret.xold = xold
    ret.x = x
    return ret

def newton(f, x0, rtol, atol, args):
    ret = NewtonResult()
    x = x0
    xold = x0*(1e0+1e30*rtol)
    fval = atol*array([1e30*ones(x0),1e30*ones(x0)])
    while(abs(fval[0]) > atol and 
          sqrt((x-xold).dot(x-xold))/(sqrt(x.dot(x))*(1e0+1e-30)) > rtol):
        fval = f(x, *args)
        xold = x
        x = x - fval[0]/fval[1]
    ret.xold = xold
    ret.x = x
    return ret
