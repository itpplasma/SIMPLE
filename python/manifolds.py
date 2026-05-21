# %%
import os
import matplotlib.pyplot as plt
import numpy as np
from common import (r0,th0,ph0,pph0,timesteps, autosave, folderforPlots,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function
from cpp_sym_midpoint import q_cpp_sym, v_cpp_sym, p_cpp_sym, dt_cpp_sym, Vtot_cpp_sym, Vpar_cpp_sym, Vper_cpp_sym, vpar_cpp_sym_vec, vper_cpp_sym_vec, En_cpp_sym, Pphi_cpp_sym, Mu_cpp_sym
from cpp_var_midpoint import q_cpp_var, v_cpp_var, p_cpp_var, dt_cpp_var, Vtot_cpp_var, Vpar_cpp_var, Vper_cpp_var, vpar_cpp_var_vec, vper_cpp_var_vec, En_cpp_var, Pphi_cpp_var, Mu_cpp_var
from cp_sym_midpoint import q_cp, v_cp, p_cp, dt_cp_sym, Vtot_cp, Vpar_cp, Vper_cp, vpar_cp_vec, vper_cp_vec, En_cp, Pphi_cp, Mu_cp
from discrete_variational import q_gc, v_gc, p_gc, dt_gc_var, Vtot_gc, Vpar_gc, Vper_gc, vpar_gc_vec, vper_gc_vec, En_gc, Pphi_gc, Mu_gc


NumInt = 4 # number of integrators to compar: CPP_sym, CPP_var, GC_var, CP_sym
f = field()
'''
def position(M, z0):

    distances = np.linalg.norm(M - z0, axis=1)

    # Find index of the minimum distance
    idx_min = np.argmin(distances)
'''
nt = timesteps




zero = 100 #th index (q(t=zero) -tupel) equal of a point on the trajectory (curve)
# Parameters of the torus
R0 = 1
z0 = 1
print(r0, th0, ph0)

#Resolution for cartesian embedding
n = 100#200 # not to be odd !!
nh = 100#100 # n/2
l = n**2 * nh # Number of Pixels



# Grid in flux coordinates

gridR = np.linspace(0,0.5,n)
gridTh = np.linspace(0,2*np.pi,n)
gridPh = np.linspace(-np.pi,np.pi,n)

a = 0.5
th, ph = np.meshgrid(gridTh, gridPh, indexing='ij')

X = (R0 + a*np.cos(th))*np.cos(ph)
Y = (R0 + a*np.cos(th))*np.sin(ph)
Z = z0 + a*np.sin(th)

# Cartesian grid
# Components in R
X = np.linspace(-2,2,n) #(0,1.8,n)
Y = np.linspace(-2,2,n) #(0,1.8,n)
Z = np.linspace(0,2,n) #(0.2,1.2,nh)

# Canonical coordinate Transformation into flux coordinates q^i
q = np.zeros([l,3])
c = 0

for x in X:
    for y in Y:
        count = -1
        R = np.sqrt((R0 - np.sqrt(x**2 + y**2))**2 + (Z - z0)**2)

        for z in Z:
            count += 1


            if R[count] <= 0.5:

                r = R[count]

                a = z-z0
                b = np.sqrt(r**2 - a**2)
                theta = np.arctan(a/b)
                cond = np.sqrt(x**2 + y**2) - R0

                if a >= 0 and cond >= 0:
                    th = theta
                elif a >= 0 and cond < 0:
                    th = np.pi - theta
                elif a < 0 and cond < 0:
                    th = np.pi - theta
                elif a < 0 and cond >= 0:
                    th = 2*np.pi + theta

                phi = np.arctan(y/x)

                if x > 0 and y > 0:
                    ph = phi
                elif x < 0 and y > 0:
                    ph = np.pi + phi
                #elif x < 0 and y < 0:
                #    ph = np.pi + phi
                elif x > 0 and y < 0:
                    ph = 2*np.pi + phi
                else: break

                q[c,:] = [r, th, ph]
                c +=1

q = q[:c,:]
print(c)

# x^i = x^i(q)
x = (R0 + q[:,0] * np.cos(q[:,1])) * np.cos(q[:,2])
y = (R0 + q[:,0] * np.cos(q[:,1])) * np.sin(q[:,2])
z = z0 + q[:,0] * np.sin(q[:,1])





#
# cp:
#
# v^x = dx/dt = dx/dq * vel^q
qcpx = []
xcpx = []
qcpy = []
xcpy = []
qcpz = []
xcpz = []
vcpx = []
vcpy = []
vcpz = []
pcpx = []
pcpy = []
pcpz = []
vparcpx = []
vparcpy = []
vparcpz = []
vpercpx = []
vpercpy = []
vpercpz = []
for i in range(len(dt_cp_sym)):

    print(i)
    traj_v_cp = v_cp[i]
    traj_q_cp = q_cp[i]
    traj_p_cp = p_cp[i]
    traj_vpar_cp = vpar_cp_vec[i]
    traj_vper_cp = vper_cp_vec[i]

    xcp_x = np.zeros([3,len(traj_q_cp[0,:])])
    vcp_x = np.zeros([3,len(traj_q_cp[0,:])])
    pcp_x = np.zeros([3,len(traj_q_cp[0,:])])
    vparcp_x = np.zeros([3,len(traj_q_cp[0,:])])
    vpercp_x = np.zeros([3,len(traj_q_cp[0,:])])

    for j in range(len(traj_q_cp[0,:])):

        A = np.array([ [np.cos(traj_q_cp[1,j])*np.cos(traj_q_cp[2,j]), -np.sin(traj_q_cp[1,j])*np.cos(traj_q_cp[2,j]), -np.sin(traj_q_cp[2,j])],
                        [np.cos(traj_q_cp[1,j])*np.sin(traj_q_cp[2,j]), -np.sin(traj_q_cp[1,j])*np.sin(traj_q_cp[2,j]),  np.cos(traj_q_cp[2,j])],
                        [               np.sin(traj_q_cp[1,j]),   (traj_q_cp[1,j]),  0] ])

        q_dot = traj_v_cp[:,j]
        qx = traj_q_cp[:,j]
        px = traj_p_cp[:,j]
        parx = traj_vpar_cp[:,j]
        perx = traj_vper_cp[:,j]

        xq = A @ qx
        x_dot = A @ q_dot
        xp = A @ px
        xpar = A @ parx
        xper = A @ perx

        xcp_x[:,j] = xq
        vcp_x[:,j] = x_dot
        pcp_x[:,j] = xp
        vparcp_x[:,j] = xpar
        vpercp_x[:,j] = xper


    # Trajectory of particle in Cartesian coordinates
    qcpx.append((R0 + traj_q_cp[0,:] * np.cos(traj_q_cp[1,:])) * np.cos(traj_q_cp[2,:]))
    qcpy.append((R0 + traj_q_cp[0,:] * np.cos(traj_q_cp[1,:])) * np.sin(traj_q_cp[2,:]))
    qcpz.append(z0 + traj_q_cp[0,:] * np.sin(traj_q_cp[1,:]))

    xcpx.append((xcp_x[0,:]))
    xcpy.append((xcp_x[1,:]))
    xcpz.append((xcp_x[2,:]))

    # Velocity in Cartesian coordinates
    vcpx.append((vcp_x[0,:])) # + qcpz
    vcpy.append((vcp_x[1,:])) # + qcpy
    vcpz.append((vcp_x[2,:])) # + qcpz

    pcpx.append((pcp_x[0,:]))
    pcpy.append((pcp_x[1,:]))
    pcpz.append((pcp_x[2,:]))

    vparcpx.append((vparcp_x[0,:]))
    vparcpy.append((vparcp_x[1,:]))
    vparcpz.append((vparcp_x[2,:]))

    vpercpx.append((vpercp_x[0,:]))
    vpercpy.append((vpercp_x[1,:]))
    vpercpz.append((vpercp_x[2,:]))

#
# gc:
#
# v^x = dx/dt = dx/dq * vel^q
qgcx = []
xgcx = []
qgcy = []
xgcy = []
qgcz = []
xgcz = []
vgcx = []
vgcy = []
vgcz = []
pgcx = []
pgcy = []
pgcz = []
vpargcx = []
vpargcy = []
vpargcz = []
vpergcx = []
vpergcy = []
vpergcz = []
for i in range(len(dt_gc_var)):

    print(i)
    traj_v_gc = v_gc[i]
    traj_q_gc = q_gc[i]
    traj_p_gc = p_gc[i]
    traj_vpar_gc = vpar_gc_vec[i]
    traj_vper_gc = vper_gc_vec[i]

    xgc_x = np.zeros([3,len(traj_q_gc[0,:])])
    vgc_x = np.zeros([3,len(traj_q_gc[0,:])])
    pgc_x = np.zeros([3,len(traj_q_gc[0,:])])
    vpar_gc_x = np.zeros([3,len(traj_q_gc[0,:])])
    vper_gc_x = np.zeros([3,len(traj_q_gc[0,:])])

    for j in range(len(traj_q_gc[0,:])):

        A = np.array([ [np.cos(traj_q_gc[1,j])*np.cos(traj_q_gc[2,j]), -np.sin(traj_q_gc[1,j])*np.cos(traj_q_gc[2,j]), -np.sin(traj_q_gc[2,j])],
                        [np.cos(traj_q_gc[1,j])*np.sin(traj_q_gc[2,j]), -np.sin(traj_q_gc[1,j])*np.sin(traj_q_gc[2,j]),  np.cos(traj_q_gc[2,j])],
                        [               np.sin(traj_q_gc[1,j]),   (traj_q_gc[1,j]),  0] ])

        q_dot = traj_v_gc[:,j]
        qx = traj_q_gc[:,j]
        px = traj_p_gc[:,j]
        parx = traj_vpar_gc[:,j]
        perx = traj_vper_gc[:,j]

        xq = A @ qx
        x_dot = A @ q_dot
        xp = A @ px
        xpar = A @ parx
        xper = A @ perx

        vgc_x[:,j] = x_dot
        xgc_x[:,j] = xq
        pgc_x[:,j] = xp
        vpar_gc_x[:,j] = xpar
        vper_gc_x[:,j] = xper


    # Trajectory of particle in Cartesian coordinates
    qgcx.append((R0 + traj_q_gc[0,:] * np.cos(traj_q_gc[1,:])) * np.cos(traj_q_gc[2,:]))
    qgcy.append((R0 + traj_q_gc[0,:] * np.cos(traj_q_gc[1,:])) * np.sin(traj_q_gc[2,:]))
    qgcz.append(z0 + traj_q_gc[0,:] * np.sin(traj_q_gc[1,:]))

    xgcx.append((xgc_x[0,:]))
    xgcy.append((xgc_x[1,:]))
    xgcz.append((xgc_x[2,:]))

    vgcx.append(vgc_x[0,:])
    vgcy.append(vgc_x[1,:])
    vgcz.append(vgc_x[2,:])

    pgcx.append(pgc_x[0,:])
    pgcy.append(pgc_x[1,:])
    pgcz.append(pgc_x[2,:])

    vpargcx.append(vpar_gc_x[0,:])
    vpargcy.append(vpar_gc_x[1,:])
    vpargcz.append(vpar_gc_x[2,:])

    vpergcx.append(vper_gc_x[0,:])
    vpergcy.append(vper_gc_x[1,:])
    vpergcz.append(vper_gc_x[2,:])

#
# cpp_symplectic:
#
# v^x = dx/dt = dx/dq * vel^q
xcppx = []
xcppy = []
xcppz = []
qcppx = []
qcppy = []
qcppz = []
vcppx = []
vcppy = []
vcppz = []
pcppx = []
pcppy = []
pcppz = []
vparcppx = []
vparcppy = []
vparcppz = []
vpercppx = []
vpercppy = []
vpercppz = []
for i in range(len(dt_cpp_sym)):

    traj_v = v_cpp_sym[i]
    traj_q = q_cpp_sym[i]
    traj_p = p_cpp_sym[i]
    traj_vpar = vpar_cpp_sym_vec[i]
    traj_vper = vper_cpp_sym_vec[i]

    xsym_x = np.zeros([3,len(traj_q[0,:])])
    vel_x = np.zeros([3,len(traj_q[0,:])])
    p_x = np.zeros([3,len(traj_q[0,:])])
    vpar_x = np.zeros([3,len(traj_q[0,:])])
    vper_x = np.zeros([3,len(traj_q[0,:])])

    for j in range(len(traj_q[0,:])):

        A = np.array([ [np.cos(traj_q[1,j])*np.cos(traj_q[2,j]), -np.sin(traj_q[1,j])*np.cos(traj_q[2,j]), -np.sin(traj_q[2,j])],
                        [np.cos(traj_q[1,j])*np.sin(traj_q[2,j]), -np.sin(traj_q[1,j])*np.sin(traj_q[2,j]),  np.cos(traj_q[2,j])],
                        [               np.sin(traj_q[1,j]),   (traj_q[1,j]),  0] ])

        q_dot = traj_v[:,j]
        qx = traj_q[:,j]
        px = traj_p[:,j]
        parx = traj_vpar[:,j]
        perx = traj_vper[:,j]

        xq = A @ qx
        x_dot = A @ q_dot
        xp = A @ px
        xpar = A @ parx
        xper = A @ perx

        xsym_x[:,j] = xq
        vel_x[:,j] = x_dot
        p_x[:,j] = xp
        vpar_x[:,j] = xpar
        vper_x[:,j] = xper

    # Trajectory of particle in Cartesian coordinates
    qcppx.append((R0 + traj_q[0,:] * np.cos(traj_q[1,:])) * np.cos(traj_q[2,:]))
    qcppy.append((R0 + traj_q[0,:] * np.cos(traj_q[1,:])) * np.sin(traj_q[2,:]))
    qcppz.append(z0 + traj_q[0,:] * np.sin(traj_q[1,:]))

    xcppx.append((xsym_x[0,:]))
    xcppy.append((xsym_x[1,:]))
    xcppz.append((xsym_x[2,:]))

    # Velocity in Cartesian coordinates
    vcppx.append(vel_x[0,:])#*50  + qcppx)
    vcppy.append(vel_x[1,:])#*50  + qcppy)
    vcppz.append(vel_x[2,:])#*50  + qcppz)

    pcppx.append(p_x[0,:])
    pcppy.append(p_x[1,:])
    pcppz.append(p_x[2,:])

    vparcppx.append(vpar_x[0,:])
    vparcppy.append(vpar_x[1,:])
    vparcppz.append(vpar_x[2,:])

    vpercppx.append(vper_x[0,:])
    vpercppy.append(vper_x[1,:])
    vpercppz.append(vper_x[2,:])


#
# cpp_variational:
#
# v^x = dx/dt = dx/dq * vel^q
xcppvarx = []
xcppvary = []
xcppvarz = []
qcppvarx = []
qcppvary = []
qcppvarz = []
vcppvarx = []
vcppvary = []
vcppvarz = []
pcppvarx = []
pcppvary = []
pcppvarz = []
vparcppvarx = []
vparcppvary = []
vparcppvarz = []
vpercppvarx = []
vpercppvary = []
vpercppvarz = []
for i in range(len(dt_cpp_var)):

    traj_v_var = v_cpp_var[i]
    traj_q_var = q_cpp_var[i]
    traj_p_var = p_cpp_var[i]
    traj_vpar_var = vpar_cpp_var_vec[i]
    traj_vper_var = vper_cpp_var_vec[i]


    xvar_x = np.zeros([3,len(traj_q_var[0,:])])
    velvar_x = np.zeros([3,len(traj_q_var[0,:])])
    pvar_x = np.zeros([3,len(traj_q_var[0,:])])
    vparvar_x = np.zeros([3,len(traj_q_var[0,:])])
    vpervar_x = np.zeros([3,len(traj_q_var[0,:])])

    for j in range(len(traj_q_var[0,:])):

        A = np.array([ [np.cos(traj_q_var[1,j])*np.cos(traj_q_var[2,j]), -np.sin(traj_q_var[1,j])*np.cos(traj_q_var[2,j]), -np.sin(traj_q_var[2,j])],
                        [np.cos(traj_q_var[1,j])*np.sin(traj_q_var[2,j]), -np.sin(traj_q_var[1,j])*np.sin(traj_q_var[2,j]),  np.cos(traj_q_var[2,j])],
                        [               np.sin(traj_q_var[1,j]),   (traj_q_var[1,j]),  0] ])

        q_dot = traj_v_var[:,j]
        qx = traj_q_var[:,j]
        px = traj_p_var[:,j]
        parx = traj_vpar_var[:,j]
        perx = traj_vper_var[:,j]

        xq = A @ qx
        x_dot = A @ q_dot
        xp = A @ px
        xpar = A @ parx
        xper = A @ perx

        xvar_x[:,j] = xq
        velvar_x[:,j] = x_dot
        pvar_x[:,j] = xp
        vparvar_x[:,j] = xpar
        vpervar_x[:,j] = xper

    # Trajectory of particle in Cartesian coordinates
    qcppvarx.append((R0 + traj_q_var[0,:] * np.cos(traj_q_var[1,:])) * np.cos(traj_q_var[2,:]))
    qcppvary.append((R0 + traj_q_var[0,:] * np.cos(traj_q_var[1,:])) * np.sin(traj_q_var[2,:]))
    qcppvarz.append(z0 + traj_q_var[0,:] * np.sin(traj_q_var[1,:]))

    xcppvarx.append((xvar_x[0,:]))
    xcppvary.append((xvar_x[1,:]))
    xcppvarz.append((xvar_x[2,:]))

    # Velocity in Cartesian coordinates
    vcppvarx.append(velvar_x[0,:])#*50  + qcppx
    vcppvary.append(velvar_x[1,:])#*50  + qcppy
    vcppvarz.append(velvar_x[2,:])#*50  + qcppz

    pcppvarx.append(pvar_x[0,:])
    pcppvary.append(pvar_x[1,:])
    pcppvarz.append(pvar_x[2,:])

    vparcppvarx.append(vparvar_x[0,:])
    vparcppvary.append(vparvar_x[1,:])
    vparcppvarz.append(vparvar_x[2,:])

    vpercppvarx.append(vpervar_x[0,:])
    vpercppvary.append(vpervar_x[1,:])
    vpercppvarz.append(vpervar_x[2,:])


'''
# Calculate Surfaces and Curves at Point p0
#
#
#
r0, th0, ph0 = traj_q[0,::100], traj_q[1,::100], traj_q[2,::100]
r0[0], th0[0], ph0[0] = traj_q[0, zero], traj_q[1, zero], traj_q[2, zero]
#
#
#
surf_R = []
surf_Th = []
surf_Phi = []
curv_R = []
curv_Th = []
curv_Phi = []

for i in range(0,len(r0)):
    print(i)
    r0i = r0[i]
    th0i = th0[i]
    ph0i = ph0[i]

    surf_Ri = np.c_[ (R0 + r0i*np.cos(q[:,1]))*np.cos(q[:,2]),
                     (R0 + r0i*np.cos(q[:,1]))*np.sin(q[:,2]),
                     z0 + r0i*np.sin(q[:,1]) ]

    surf_Thi = np.c_[ (R0 + q[:,0]*np.cos(th0i))*np.cos(q[:,2]),
                      (R0 + q[:,0]*np.cos(th0i))*np.sin(q[:,2]),
                      z0 + q[:,0]*np.sin(th0i) ]

    surf_Phii = np.c_[ (R0 + q[:,0]*np.cos(q[:,1]))*np.cos(ph0i),
                       (R0 + q[:,0]*np.cos(q[:,1]))*np.sin(ph0i),
                       z0 + q[:,0]*np.sin(q[:,1]) ]

    curv_Ri = np.c_[ (R0 + q[:,0]*np.cos(th0i))*np.cos(ph0i),
                      (R0 + q[:,0]*np.cos(th0i))*np.sin(ph0i),
                      z0 + q[:,0]*np.sin(th0i) ]

    curv_Thi = np.c_[ (R0 + r0i*np.cos(q[:,1]))*np.cos(ph0i),
                       (R0 + r0i*np.cos(q[:,1]))*np.sin(ph0i),
                       z0 + r0i*np.sin(q[:,1]) ]

    curv_Phii = np.c_[ (R0 + r0i*np.cos(th0i))*np.cos(q[:,2]),
                       (R0 + r0i*np.cos(th0i))*np.sin(q[:,2]),
                       z0 + r0i*np.sin(th0i)*np.ones(np.size(q[:,2])) ]

    surf_R.append(surf_Ri)
    surf_Th.append(surf_Thi)
    surf_Phi.append(surf_Phii)
    curv_R.append(curv_Ri)
    curv_Th.append(curv_Thi)
    curv_Phi.append(curv_Phii)
'''

# Interpolation of the field on the grid points
# or
# Field evaluation on trajectory at grid points
#
# Input: q --- random grid of (possible) configurations in flux coordinates
#        zero --- index of point of interest
#        traj_q --- trajectory of particle
#        c --- number of points q

def static_B_field(q, zero, traj_q, c):

    B = np.zeros(c)
    k = 0

    for i in range(0,c):
        r = q[i,0]
        th = q[i,1]
        ph = q[i,2]
        f.evaluate(r,th,ph)
        B[i] = f.B


        if i == 0:
            const_B = []
            f.evaluate(traj_q[0,zero], traj_q[1,zero], traj_q[2,zero])
            const_B.append([traj_q[0,zero], traj_q[1,zero], traj_q[2,zero], f.B])
            eps = abs(const_B[i][3] - B[0])
            print(eps)

        elif eps <= 1e-3:
            const_B.append([q[i-1,0], q[i-1,1], q[i-1,2], B[i-1]])
            eps = abs(B[i] - const_B[0][3])
            k += 1

        else:
            eps = abs(const_B[0][3] - B[i])



    print('c =' , c)
    print('k =' , k)
    # with k the number of points with almost the same B field value

    arr = np.array(const_B)
    rB = arr[:,0]
    thB = arr[:,1]
    phB = arr[:,2]
    B_const = arr[:,3]
    print(len(rB))

    consBx = (R0 + rB * np.cos(thB)) * np.cos(phB)
    consBy = (R0 + rB * np.cos(thB)) * np.sin(phB)
    consBz = z0 + rB * np.sin(thB)
    consB = np.array([consBx, consBy, consBz])

    # Field evaluation on trajectory
    B_traj = np.zeros(np.size(traj_q[0,:]))
    for j in range(0,np.size(traj_q[0,:])):
        r = traj_q[0,j]
        th = traj_q[1,j]
        ph = traj_q[2,j]
        f.evaluate(r,th,ph)
        B_traj[j] = f.B


    return B, consB, B_const, B_traj

# Output: B --- field values on grid points q
#         consB --- cartesian coordinates of constant field for point of interest
#         B_const --- flux coordinates -..-
#         B_traj --- field values on traj.

#B, consB, B_const, B_traj = static_B_field(q, zero, traj_q, c)






'''
#############################################################


#plot magnetic field axis lines

metric = lambda z: {
    "_ii":  np.array([1, z[0]**2, (1 + z[0]*np.cos(z[1]))**2]),
    "^ii":  np.array([1, 1/z[0]**2, 1/(1 + z[0]*np.cos(z[1]))**2]),
    "d_11": np.array([0,0,0]),
    "d_22": np.array([2*z[0], 0,0]),
    "d_33": np.array([2*(1+z[0]*np.cos(z[1]))*np.cos(z[1]), -2*(1+z[0]*np.cos(z[1]))*np.sin(z[1]), 0]),
}

noc = 100 # number of curves on traj.
hath = np.zeros([3,noc])
hath0 = np.zeros([3,noc])

R0 = 1
z0 = 1



for j in range(noc):

    j100 = j*10


    hath0[0,j] = (R0 + traj_q_cp[0,j100] * np.cos(traj_q_cp[1,j100])) * np.cos(traj_q_cp[2,j100])
    hath0[1,j] = (R0 + traj_q_cp[0,j100] * np.cos(traj_q_cp[1,j100])) * np.sin(traj_q_cp[2,j100])
    hath0[2,j] = z0 + traj_q_cp[0,j100] * np.sin(traj_q_cp[1,j100])

    f.evaluate(traj_q_cp[0,j100], traj_q_cp[1,j100], traj_q_cp[2,j100])
    g = metric(traj_q_cp[:,j100])

    co_h = np.array([f.co_hr, f.co_hth, f.co_hph])
    con_h = g['^ii']*co_h

    hath[0,j] = con_h[0] * ( np.cos(traj_q_cp[1,j100]) * np.cos(traj_q_cp[2,j100]) ) - con_h[1] * (np.sin(traj_q_cp[1,j100]) * np.cos(traj_q_cp[2,j100]) ) - con_h[2] * ( np.sin(traj_q_cp[2,j100]) )
    hath[1,j] = con_h[0] * ( np.cos(traj_q_cp[1,j100]) * np.sin(traj_q_cp[2,j100]) ) - con_h[1] * (np.sin(traj_q_cp[1,j100]) * np.sin(traj_q_cp[2,j100]) ) + con_h[2] * ( np.cos(traj_q_cp[2,j100]) )
    hath[2,j] = con_h[0] * ( np.sin(traj_q_cp[1,j100]) ) + con_h[1] * (np.cos(traj_q_cp[1,j100]) )

    #hath[0,j] = (R0 + con_h[0] * np.cos(con_h[1])) * np.cos(con_h[2])
    #hath[1,j] = (R0 + con_h[0] * np.cos(con_h[1])) * np.sin(con_h[2])
    #hath[2,j] = z0 + con_h[0] * np.sin(con_h[1])
    #hath[0,j] = f.co_hr*hath0[0,j]
    #hath[1,j] = f.co_hth*hath0[1,j]
    #hath[2,j] = f.co_hph*hath0[2,j]



print('hath0 :', hath0[:,0])
print('hath :', hath[:,0])
print(hath)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(qcpx[:3],qcpy[:3],qcpz[:3], label='r cp')
ax.scatter(qcpx, qcpy, qcpz, c="r", s=1)

#ax.plot([hath0[0,0], hath[0,0]], [hath0[1,0], hath[1,0]], [hath0[2,0], hath[2,0]], color='green', label='Magnetic Field Line')
ax.scatter(hath0[0,:], hath0[1,:], hath0[2,:], color='blue', marker='x', label='Magnetic Field Line Start')

ax.scatter(hath[0,:]*0.1+hath0[0,:], hath[1,:]*0.1+hath0[1,:], hath[2,:]*0.1+hath0[2,:], color='green', marker='o', label='Magnetic Field Line End')
for j in range(noc):
    ax.plot([hath0[0,j], 0.1*hath[0,j]+hath0[0,j]], [hath0[1,j], 0.1*hath[1,j]+hath0[1,j]], [hath0[2,j], 0.1*hath[2,j]+hath0[2,j]], color='green', label='Magnetic Field Line')


ax.plot(qcppx,qcppy,qcppz, c='g', label='1 cpp')
#ax.plot(qcppx[10:20],qcppy[10:20],qcppz[10:20], c='b', label='2 cpp')
#ax.plot(qcppx[20:30],qcppy[20:30],qcppz[20:30], c='orange', label='3 cpp')

ax.scatter(vcppx[::10]*80+qcppx[::10], vcppy[::10]*80+qcppy[::10], vcppz[::10]*80+qcppz[::10], c='k', s=0.1, label='cpp velocity')

for j in range(noc):
    ax.plot([qcppx[j*10], vcppx[j*10]*80+qcppx[j*10]], [qcppy[j*10], vcppy[j*10]*80+qcppy[j*10]], [qcppz[j*10], vcppz[j*10]*80+qcppz[j*10]], c='y')



for j in range(noc):
    ax.plot([qcpx[j*10], 10*vcpx[j*10]+qcpx[j*10]], [qcpy[j*10], 10*vcpy[j*10]+qcpy[j*10]], [qcpz[j*10], 10*vcpz[j*10]+qcpz[j*10]], c='b')


plt.show()
'''


#%%
# Plotting of energy error for different dts. sim_time : nt*dt = const.
# -> nt = sim_time/dt (-add to each integrator)
# In common we use:
#dt_cp_sym = [0.01, 0.1, 1, 10]
#dt_cpp_sym = [0.01, 0.1, 1, 10]
#dt_cpp_var = [0.01, 0.1, 1, 10]
#dt_gc_var = [0.01, 0.1, 1, 10]
#timesteps = 100 # = sim_time in that case

E = [En_cp, En_cpp_sym, En_cpp_var, En_gc]
DT = [dt_cp_sym, dt_cpp_sym, dt_cpp_var, dt_gc_var]

eerr_cp = []
eerr_cpp_sym = []
eerr_cpp_var = []
eerr_gc_var = []

for i in range(len(dt_cp_sym)):
    rel_cp = (En_cp[i] - En_cp[i][0]) / En_cp[i][0]
    rel_cpp_sym = (En_cpp_sym[i] - En_cpp_sym[i][0]) / En_cpp_sym[i][0]
    rel_cpp_var = (En_cpp_var[i] - En_cpp_var[i][0]) / En_cpp_var[i][0]
    rel_gc_var = (En_gc[i] - En_gc[i][0]) / En_gc[i][0]

    eerr_cp.append(np.max(np.abs(rel_cp)))
    eerr_cpp_sym.append(np.max(np.abs(rel_cpp_sym)))
    eerr_cpp_var.append(np.max(np.abs(rel_cpp_var)))
    eerr_gc_var.append(np.max(np.abs(rel_gc_var)))

eerr_cp = np.array(eerr_cp)
eerr_cpp_sym = np.array(eerr_cpp_sym)
eerr_cpp_var = np.array(eerr_cpp_var)
eerr_gc_var = np.array(eerr_gc_var)


plt.figure(figsize=(7, 5))

plt.loglog(np.array(dt_cp_sym), np.array(eerr_cp), 'o-', label='CP sym')
plt.loglog(np.array(dt_cpp_sym), np.array(eerr_cpp_sym), 's-', label='CPP sym')
plt.loglog(np.array(dt_cpp_var), np.array(eerr_cpp_var), 'd-', label='CPP var')
plt.loglog(np.array(dt_gc_var), np.array(eerr_gc_var), 'x-', label='GC var')


# reference line ~ dt^2
dt_ref = np.array(dt_cpp_sym)
err_ref = eerr_cpp_sym[-1] * (dt_ref / dt_ref[-1])**2
plt.loglog(dt_ref, err_ref, 'k--', label=r'$\sim \Delta t^2$')

plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\max_t |\Delta E/E_0|$')
plt.grid(True, which='both')
plt.legend()

if autosave == 1:
    plt.savefig(os.path.join(folderforPlots, 'energy_error.png'), dpi=300)

plt.show()




#%%
#from matplotlib.ticker import MaxNLocator
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

### PLOTTING ###

#  Velocity evolution
Int = ['cp\_sym', 'cpp\_sym', 'cpp\_var', 'gc\_var']
VTOT = [Vtot_cp, Vtot_cpp_sym, Vtot_cpp_var, Vtot_gc]
VPAR = [Vpar_cp, Vpar_cpp_sym, Vpar_cpp_var, Vpar_gc]
VPER = [Vper_cp, Vper_cpp_sym, Vper_cpp_var, Vper_gc]
DT = [dt_cp_sym, dt_cpp_sym, dt_cpp_var, dt_gc_var]

k = nt


for j in range(NumInt):
    vtot = VTOT[j]
    vpar = VPAR[j]
    vper = VPER[j]
    dt = DT[j]

    fig = plt.figure(figsize=(16,5))

    for i in range(len(dt)):
        #k = int(nt/dt[i]) #//100
        ax = fig.add_subplot(1, len(dt), i+1)
        t = np.linspace(0, dt[i]*nt, nt+1)

        ax.scatter(t[:k], vpar[i][:k], s=1, marker='o', c='steelblue')
        ax.scatter(t[:k], vper[i][:k], s=1, marker='x', c='darkorange')

        ax.plot(t[:k], vtot[i][:k], lw=1, c='seagreen', label=r'$v$')
        ax.plot(t[:k], vpar[i][:k], lw=0.2, c='steelblue', label=r'$v_\parallel$')
        ax.plot(t[:k], vper[i][:k], lw=0.2, c='darkorange', label=r'$v_\perp$')

        ax.set_title(r'$\Delta t_{\rm ' + str(Int[j]) + r'} = ' + str(dt[i]) + r'$')

        xmax = min(40000, t[-1])
        ax.set_xlim(0, xmax)
        ax.set_xticks([0, xmax/2, xmax])
        ax.set_yticks(ax.get_yticks()[::2])

        ax.set_xlabel(r'$t$')

        if i == 0:
            ax.set_ylabel(r'$v$')
            leg = ax.legend(loc='center right', markerscale=2, fontsize=12)
            for line in leg.get_lines():
                line.set_linewidth(2)

    if autosave == 1:
        fig.savefig(os.path.join(folderforPlots, 'vel_dt ' + str(Int[j]) + '.png'), dpi=300)

    plt.show()

#%%
# ??????????????????????????
#for i in range(NumInt):
#    dt = DT[i]
#
#    for j in range(len(dt)):
#        fig = plt.figure(figsize=(16,5))
#        ax = fig.add_subplot(1, len(dt), j+1)
#        t = np.linspace(0, dt[j]*nt, nt+1)
#
#        ax.scatter(t[:k], B_traj[i][:k], s=1, marker='o', c='steelblue')
#        ax.set_title(r'$\Delta t_{\rm ' + str(Int[i]) + r'} = ' + str(dt[j]) + r'$')
#
#        xmax = min(40000, t[-1])
#        ax.set_xlim(0, xmax)
#        ax.set_xticks([0, xmax/2, xmax])
#        ax.set_yticks(ax.get_yticks()[::2])
#
#        ax.set_xlabel(r'$t$')
#        ax.set_ylabel(r'$B$')

#%%
E = [En_cp, En_cpp_sym, En_cpp_var, En_gc]
PPhi = [Pphi_cp, Pphi_cpp_sym, Pphi_cpp_var, Pphi_gc]
MU = [Mu_cp, Mu_cpp_sym, Mu_cpp_var, Mu_gc]

row_labels = [
    r'$\frac{\Delta E}{E_0}$',
    r'$\frac{\Delta \mu}{\mu_0}$',
    r'$\frac{\Delta p_\phi}{p_{\phi,0}}$'
]

for i in range(NumInt):
    En = E[i]
    mui = MU[i]
    pphii = PPhi[i]
    dt = DT[i]

    fig = plt.figure(figsize=(16,12))

    for j in range(len(dt)):
        #k = int(nt/dt[i])

        t = np.linspace(0, dt[j]*nt, nt+1)
        ax = fig.add_subplot(3,  len(dt), j+1)
        ax.scatter(t[:k], (En[j][:k]-En[j][0])/(En[j][0]), s=1, marker='o', c='firebrick', label=r'$\frac{\Delta E}{E_0}$')
        if j == 0:
            ax.set_ylabel(row_labels[0], fontsize=16)

        ax = fig.add_subplot(3,  len(dt), j+1+len(dt))
        ax.scatter(t[:k], (mui[j][:k]-1e-5)/(1e-5), s=1, marker='x', c='black', label=r'$\frac{\Delta \mu}{\mu_0}$')
        if j == 0:
            ax.set_ylabel(row_labels[1], fontsize=16)

        ax = fig.add_subplot(3,  len(dt), j+1+2*len(dt))
        ax.scatter(t[:k], (pphii[j][:k]-pphii[j][0])/(pphii[j][0]), s=1, marker='^', c='seagreen', label=r'$\frac{\Delta p_\phi}{p_{\phi,0}}$')
        ax.set_xlabel(r'$\Delta t = ' + str(dt[j]) + r'$')
        if j == 0:
            ax.set_ylabel(row_labels[2], fontsize=16)

    fig.suptitle(r'Conservation of Invariants for ${\rm ' + str(Int[i]) + r'}$', fontsize=20, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if autosave == 1:
        fig.savefig(os.path.join(folderforPlots, 'invariants_dt ' + str(Int[i]) + '.png'), dpi=300)


plt.show()

#%%

t = [1, ]

#%%

# plot Banana orbit (Curie plot idk)

fig, ax = plt.subplots()

plot_orbit(q_cpp_sym[0], ax=ax)
plot_orbit(q_cp[0], ax=ax)
#plot_orbit(consB, ax=ax)
plt.show()

fig = plt.figure(figsize=(20, 6))

# 1) Charged particle
ax = fig.add_subplot(1, 4, 1)
ax.scatter(
    q_cp[-1][0, :] * np.cos(q_cp[-1][1, :]),
    q_cp[-1][0, :] * np.sin(q_cp[-1][1, :]),
    s=0.1,
    c='r',
    label='Charged Particle'
)
ax.set_title('Charged Particle')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(markerscale=20, fontsize=10)

# 2) CPP symplectic
ax = fig.add_subplot(1, 4, 2)
ax.scatter(
    q_cpp_sym[-1][0, :] * np.cos(q_cpp_sym[-1][1, :]),
    q_cpp_sym[-1][0, :] * np.sin(q_cpp_sym[-1][1, :]),
    s=0.1,
    c='b',
    label='CPP - Symplectic'
)
ax.set_title('CPP - Symplectic')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(markerscale=20, fontsize=10)

# 3) CPP variational
ax = fig.add_subplot(1, 4, 3)
ax.scatter(
    q_cpp_var[-1][0, :] * np.cos(q_cpp_var[-1][1, :]),
    q_cpp_var[-1][0, :] * np.sin(q_cpp_var[-1][1, :]),
    s=0.1,
    c='g',
    label='CPP - Variational'
)
ax.set_title('CPP - Variational')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(markerscale=20, fontsize=10)

# 4) Guiding center
ax = fig.add_subplot(1, 4, 4)
ax.scatter(
    q_gc[-1][0, :] * np.cos(q_gc[-1][1, :]),
    q_gc[-1][0, :] * np.sin(q_gc[-1][1, :]),
    s=0.1,
    c='purple',
    label='Guiding Center'
)
ax.set_title('Guiding Center')
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(markerscale=20, fontsize=10)

fig.subplots_adjust(
    left=0.04,
    right=0.98,
    bottom=0.12,
    top=0.88,
    wspace=0.25
)

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'comparison_banana.png'), dpi=300)

plt.show()

#%%

fig = plt.figure(figsize=(20, 6))

for j in range(len(dt_gc_var)):
    ax = fig.add_subplot(1, len(dt_gc_var), j+1)
    ax.scatter(q_cpp_sym[j][0,:]*np.cos(q_cpp_sym[j][1,:]),
               q_cpp_sym[j][0,:]*np.sin(q_cpp_sym[j][1,:]), s=0.1, c='b', label='CPP - Symplectic')
    ax.scatter(q_cpp_var[j][0,:]*np.cos(q_cpp_var[j][1,:]),
               q_cpp_var[j][0,:]*np.sin(q_cpp_var[j][1,:]), s=0.1, c='g', label='CPP - Variational')
    ax.plot(q_cp[2][0,:]*np.cos(q_cp[2][1,:]),
               q_cp[2][0,:]*np.sin(q_cp[2][1,:]), lw=0.02, c='r', label='Charged Particle for dt = ' + str(dt_cp_sym[2]))
    ax.scatter(q_gc[j][0,:]*np.cos(q_gc[j][1,:]),
               q_gc[j][0,:]*np.sin(q_gc[j][1,:]), s=0.1, marker='o', c='purple', label='Guiding Center Particle')

    plt.title(r"$\Delta t = " + str(dt_gc_var[j]) + "$")
    if j == 0:
        leg = ax.legend(loc='center left', markerscale=20, fontsize=10)
        for line in leg.get_lines():
            line.set_linewidth(2)

if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'banana.png'), dpi=300)
plt.show()


#%%

# embedding visualization in 3D euclidean space

#%%
# Config. space and Tangent bundle
def setup_3d_axis(
    ax,
    xlim,
    ylim,
    zlim,
    elev,
    azim,
    axis_off=False,
):
    #ax.set_xlim(*xlim)
    #ax.set_ylim(*ylim)
   # ax.set_zlim(*zlim)

    ax.view_init(elev=elev, azim=azim)
    ax.grid(False)

    if axis_off:
        ax.set_axis_off()
    else:
        #ax.set_xticks([])
        ax.set_xlabel('x')
        ax.set_xticks(ax.get_xticks()[::2])
        ax.set_ylabel('y')
        ax.set_yticks(ax.get_yticks()[::2])
        ax.set_zticks([])



def plot_tangent_bundle_comparison(
    qx, qy, qz,
    vx, vy, vz,
    px, py, pz,
    dts,
    label,
    qcpx_ref, qcpy_ref, qcpz_ref,
    autosave=0,
    folderforPlots=None,
    filename=None,
    v_amp=1,#(2*np.pi)**2,
    quiver_step=12,
    figsize=None,
):


    n_col = len(dts)

    if figsize is None:
        figsize = (4 * n_col, 10)

    fig = plt.figure(figsize=figsize)

    # Reference limits from charged particle trajectory
    xlim_ref = (qcpx_ref.min(), qcpx_ref.max())
    ylim_ref = (qcpy_ref.min(), qcpy_ref.max())
    zlim_ref = (qcpz_ref.min(), qcpz_ref.max())

    # Limits from current model
    xlim_model = (min([arr.min() for arr in qx]), max([arr.max() for arr in qx]))
    ylim_model = (min([arr.min() for arr in qy]), max([arr.max() for arr in qy]))
    zlim_model = (min([arr.min() for arr in qz]), max([arr.max() for arr in qz]))


    xlim_model = (0, 2.5)


    for j in range(n_col):




        qX = qx[j]
        qY = qy[j]
        qZ = qz[j]

        vX = vx[j]
        vY = vy[j]
        vZ = vz[j]

        pX = px[j]
        pY = py[j]
        pZ = pz[j]


        v_amp = dts[j]*np.sqrt((2*np.pi))
        p_amp = np.sqrt((2*np.pi)**3/3) #dts[0]/dts[j]#(2*np.pi)
        print(r)

        # ------------------------------------------------------------
        # Row 1: Configuration-space trajectory
        # ------------------------------------------------------------
        ax = fig.add_subplot(3, n_col, j + 1 , projection="3d")

        #ax.plot(qcpx_ref, qcpy_ref, qcpz_ref, c="r", lw=0.15)
        #ax.plot(qX, qY, qZ, c="b", lw=0.01)
        ax.scatter(qX, qY, qZ, c="b", s=0.02)

        setup_3d_axis(
            ax,
            xlim_model,
            ylim_model,
            zlim_model,
            elev=-90,
            azim=0
        )

        # ------------------------------------------------------------
        # Row 2: Tangent bundle visualization
        # ------------------------------------------------------------
        ax = fig.add_subplot(3, n_col, n_col + j + 1, projection="3d")

        #ax.plot(qcpx_ref, qcpy_ref, qcpz_ref, c="r", lw=0.15)
        ax.plot(qX, qY, qZ, c="b", lw=0.2)
        #ax.scatter(qX, qY, qZ, c="b", s=0.02)

        ax.scatter(
            qX + v_amp * vX,
            qY + v_amp * vY,
            qZ + v_amp * vZ,
            c="g",
            s=0.02
        )
        #ax.plot(qX+ v_amp * vX, qY + v_amp * vY, qZ + v_amp * vZ, c="g", lw=0.01)

        #ax.quiver(
        #    qX[::quiver_step],
        #    qY[::quiver_step],
        #    qZ[::quiver_step],
        #    vX[::quiver_step],
        #    vY[::quiver_step],
        #    vZ[::quiver_step],
        #    length=v_amp,
        #    normalize=False,
        #    color="g",
        #    lw=0.6
        #)

        setup_3d_axis(
            ax,
            xlim_model,
            ylim_model,
            zlim_model,
            elev=-90,
            azim=0
        )

        # ------------------------------------------------------------
        # Row 3: Alternative / zoomed view
        # ------------------------------------------------------------
        ax = fig.add_subplot(3, n_col, j + 1 + 2 * n_col, projection="3d")

        #ax.plot(qcpx_ref, qcpy_ref, qcpz_ref, c="r", lw=0.15)
        #ax.scatter(qX, qY, qZ, c="b", s=0.02)
        ax.plot(qX, qY, qZ, c="b", lw=0.2)

        ax.scatter(
            qX + p_amp * pX,
            qY + p_amp * pY,
            qZ + p_amp * pZ,
            c="r",
            s=0.02
        )

        setup_3d_axis(
            ax,
            xlim_model,
            ylim_model,
            zlim_model,
            elev=-90,
            azim=0
        )

        # Column labels
        ax.text2D(
            0.5,
            0.05,
            rf"$\Delta t_{{\rm {label}}} = {dts[j]}$",
            transform=ax.transAxes,
            ha="center",
            fontsize=12
        )

    fig.text(
        0.5,
        0.6,
        rf"$R(q) + \sqrt{{2\pi}} \cdot \Delta t \cdot \, v(q), \qquad $ amp = {v_amp}",
        ha="center",
        fontsize=14
    )

    fig.text(
        0.5,
        0.96,
        r"Trajectory in configuration space $Q$",
        ha="center",
        fontsize=14
    )

    fig.text(
        0.5,
        0.25,
        rf"$R(q) + \sqrt{{(2\pi)^3}} \cdot \, p(q), \qquad $ amp = {p_amp}",
        ha="center",
        fontsize=14
    )

    plt.subplots_adjust(
        left=0.02,
        right=0.98,
        top=0.98,
        bottom=-0.1,
        wspace=-0.1,
        hspace=0.3
    )

    if autosave == 1 and folderforPlots is not None:
        if filename is None:
            filename = f"{label}_xq.png"

        fig.savefig(
            os.path.join(folderforPlots, filename),
            dpi=300,
            bbox_inches="tight"
        )

    plt.show()


deltt = 0

qcpx0 = qcpx[deltt]
qcpy0 = qcpy[deltt]
qcpz0 = qcpz[deltt]



plot_tangent_bundle_comparison(
    qgcx, qgcy, qgcz,
    vgcx, vgcy, vgcz,
    pgcx, pgcy, pgcz,
    dt_gc_var,
    label="gc\_var",
    qcpx_ref=qcpx0,
    qcpy_ref=qcpy0,
    qcpz_ref=qcpz0,
    autosave=autosave,
    folderforPlots=folderforPlots,
    filename="gc_var_xq.png"
)



plot_tangent_bundle_comparison(
    qcpx, qcpy, qcpz,
    vcpx, vcpy, vcpz,
    pcpx, pcpy, pcpz,
    dt_cp_sym,
    label="cp\_sym",
    qcpx_ref=qcpx0,
    qcpy_ref=qcpy0,
    qcpz_ref=qcpz0,
    autosave=autosave,
    folderforPlots=folderforPlots,
    filename="cp_sym_xq.png"
)


plot_tangent_bundle_comparison(
    qcppx, qcppy, qcppz,
    vcppx, vcppy, vcppz,
    pcppx, pcppy, pcppz,
    dt_cpp_sym,
    label="cpp\_sym",
    qcpx_ref=qcpx0,
    qcpy_ref=qcpy0,
    qcpz_ref=qcpz0,
    autosave=autosave,
    folderforPlots=folderforPlots,
    filename="cpp_sym_xq.png"
)

plot_tangent_bundle_comparison(
    qcppvarx, qcppvary, qcppvarz,
    vcppvarx, vcppvary, vcppvarz,
    pcppvarx, pcppvary, pcppvarz,
    dt_cpp_var,
    label="cpp\_var",
    qcpx_ref=qcpx0,
    qcpy_ref=qcpy0,
    qcpz_ref=qcpz0,
    autosave=autosave,
    folderforPlots=folderforPlots,
    filename="cpp_var_xq.png"
)

#%%
'''

fig = plt.figure(figsize=(18,12))
n_col = len(dt_cpp_sym)


for j in range(len(dt_cpp_sym)):
    ax1 = fig.add_subplot(len(dt_cp_sym), n_col, 3*n_col+j+1, projection='3d')
    ax1.plot(xcpx[j], xcpy[j], xcpz[j], color='k', lw=0.01)
    ax1.scatter(xcpx[j], xcpy[j], xcpz[j], color='gray', s=0.01)


    slice_ids = np.linspace(0, nt, nt, dtype=int)
    all_idx = []

    for k in slice_ids:
        if k < len(xcppy[j]):
            y0 = xcppy[j][k]
            idx = np.where(np.abs(xcppy[j] - y0) < 1e-3)[0]
            all_idx.extend(idx)
            ax1.plot(xcpx[j][idx], xcpy[j][idx], xcpz[j][idx], color='b', lw=0.01)
            ax1.scatter(xcpx[j][idx], xcpy[j][idx], xcpz[j][idx], color='violet', s=0.05, marker='o')
    #ax1.set_title(r'$\Delta t_{{\rm cp}} = ' + str(dt_cp_sym[j]) + '$')
    ax1.view_init(elev=-90, azim=0)
    #dxj = 0.5 * dt_cp_sym[j]
    #ax1.set_ylim(-dxj,dxj)
    #ax1.set_box_aspect([4, 1, 1])
    if dt_cp_sym[j] != 10:
        ax1.set_axis_off()

    ax2 = fig.add_subplot(len(dt_cpp_sym), n_col, n_col+1+j , projection='3d')
    ax2.set_facecolor('none')
    #ax2.plot(xcppx[j], xcppy[j], xcppz[j], color='k', lw=0.01)
    ax2.scatter(xcppx[j], xcppy[j], xcppz[j], color='gray', s=0.01)

    slice_ids = np.linspace(0, nt, nt, dtype=int)
    all_idx = []

    for k in slice_ids:
        if k < len(xcppy[j]):
            y0 = xcppy[j][k]
            idx = np.where(np.abs(xcppy[j] - y0) < 1e-3)[0]
            all_idx.extend(idx)
            ax2.plot(xcppx[j][idx], xcppy[j][idx], xcppz[j][idx], color='b', lw=0.01)
            ax2.scatter(xcppx[j][idx], xcppy[j][idx], xcppz[j][idx], color='violet', s=0.05, marker='o')

    #all_idx = np.unique(all_idx)

    ax2.view_init(elev=90, azim=0)
    dxj = 0.5 * dt_cpp_sym[j]
    ax2.set_ylim(-dxj, dxj)
    if dt_cpp_sym[j] != 10:
        ax2.set_axis_off()



    ax3 = fig.add_subplot(len(dt_cpp_sym), n_col, 2*n_col+1+j, projection='3d')
    ax3.set_facecolor('none')
    #ax3.plot(xcppvarx[j], xcppvary[j], xcppvarz[j], color='k', lw=0.01)
    ax3.scatter(xcppvarx[j], xcppvary[j], xcppvarz[j], color='gray', s=0.01)
    print('nt = ', nt)
    slice_ids = np.linspace(0, nt, nt, dtype=int)
    all_idx = []

    for k in slice_ids:
        if k < len(xcppvary[j]):
            y0 = xcppvary[j][k]
            idx = np.where(np.abs(xcppvary[j] - y0) < 1e-3)[0]
            all_idx.extend(idx)
            ax3.plot(xcppvarx[j][idx], xcppvary[j][idx], xcppvarz[j][idx], color='b', lw=0.01)
            ax3.scatter(xcppvarx[j][idx], xcppvary[j][idx], xcppvarz[j][idx], color='violet', s=0.05, marker='o')

    #all_idx = np.unique(all_idx)

    ax3.view_init(elev=0, azim=0)
    dxj = 0.5 * dt_cpp_var[j]
    ax3.set_ylim(-dxj, dxj)
    ax3.set_axis_off()
    ax3.axhline(0.1, color='k', lw=5)
    #
    ## inset histogram inside ax3
    #
    #axins = inset_axes(ax3, width="100%", height="100%", loc="center")
    ## make inset axis transparent
    #axins.set_facecolor('none')
    #axins.patch.set_alpha(0)
    ## remove frame
    #for spine in axins.spines.values():
    #    spine.set_visible(False)

    #axins.hist(xcppvarz[j][all_idx], bins=200, alpha=0.3, color='b', edgecolor='b')
    #axins.set_xticks([])
    #axins.set_yticks([])



    ax4 = fig.add_subplot(len(dt_gc_var), n_col, 1+j, projection='3d')
    ax4.set_facecolor('none')
    #ax4.plot(xgcx[j], xgcy[j], xgcz[j], color='k', lw=0.01)
    ax4.scatter(xgcx[j], xgcy[j], xgcz[j], color='gray', s=0.01)

    slice_ids = np.linspace(0, nt, nt, dtype=int)
    all_idx = []

    for k in slice_ids:
        if k < len(xgcy[j]):
            x0 = xgcx[j][k]
            idx = np.where(np.abs(xgcx[j] - x0) < 1e-3)[0]
            all_idx.extend(idx)
            ax4.plot(xgcx[j][idx], xgcy[j][idx], xgcz[j][idx], color='b', lw=0.01)
            ax4.scatter(xgcx[j][idx], xgcy[j][idx], xgcz[j][idx], color='violet', s=0.05, marker='o')

    #all_idx = np.unique(all_idx)
    ax4.view_init(elev=0, azim=0)
    dxj = 0.5 * dt_cpp_var[j]
    ax4.set_ylim(-dxj, dxj)
    ax4.set_axis_off()



col_titles_top = [rf'$\Delta t_{{\rm cpp}}$ = {dt_cpp_sym[0]}',rf'$\Delta t_{{\rm cp}}$ = {dt_cpp_sym[1]}', rf'$\Delta t_{{\rm cpp}}$ = {dt_cpp_sym[2]}', rf'$\Delta t_{{\rm cpp}}$ = {dt_cpp_sym[3]}']
col_titles_bot = [rf'$\Delta t_{{\rm cp}}$ = {dt_cp_sym[0]}',rf'$\Delta t_{{\rm cp}}$ = {dt_cp_sym[1]}', rf'$\Delta t_{{\rm cp}}$ = {dt_cp_sym[2]}', rf'$\Delta t_{{\rm cp}}$ = {dt_cp_sym[3]}']
xcols = [0.16, 0.38, 0.62, 0.84]

for x, title in zip(xcols, col_titles_top):
    fig.text(x, 0.985, title, ha='center', va='top', fontsize=16)
for x, title in zip(xcols, col_titles_bot):
    fig.text(x, 0.01, title, ha='center', va='bottom', fontsize=16)

# row labels on left and right
row_titles =['Guiding Center', 'CPP Symplectic', 'CPP Variational', 'Classical Particle']
for i, title in enumerate(row_titles):
    y = 1 - (i + 0.5) / len(dt_cpp_sym)

    fig.text(0.98, y,
            title,
            ha='right', va='center', fontsize=12)

plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.05, wspace=-0.2, hspace=-0.3)
if autosave == 1:
    fig.savefig(os.path.join(folderforPlots, 'rel_position_x(q(x)).png'), dpi=300)
plt.show()

'''

#
#
#
###
#
#
#
###
#
#
#

#%%
from matplotlib.lines import Line2D



def turning_point_distances(x, y, z, vx, vy, vz):
    """
    Finds distances between consecutive points where vx, vy, vz change sign.

    Returns:
        dx_turn: x-distance between points where vx changes sign
        dy_turn: y-distance between points where vy changes sign
        dz_turn: z-distance between points where vz changes sign

        d3_vx: 3D distance between points where vx changes sign
        d3_vy: 3D distance between points where vy changes sign
        d3_vz: 3D distance between points where vz changes sign
    """

    def sign_change_indices(v):
        s = np.sign(v).copy()

        for i in range(1, len(s)):
            if s[i] == 0:
                s[i] = s[i-1]

        return np.where(s[:-1] * s[1:] < 0)[0]

    idx_vx = sign_change_indices(vx)
    idx_vy = sign_change_indices(vy)
    idx_vz = sign_change_indices(vz)

    # points where velocity components change sign
    px_vx = x[idx_vx]
    py_vx = y[idx_vx]
    pz_vx = z[idx_vx]

    px_vy = x[idx_vy]
    py_vy = y[idx_vy]
    pz_vy = z[idx_vy]

    px_vz = x[idx_vz]
    py_vz = y[idx_vz]
    pz_vz = z[idx_vz]

    # coordinate distances
    dx_turn = np.abs(np.diff(px_vx))
    dy_turn = np.abs(np.diff(py_vy))
    dz_turn = np.abs(np.diff(pz_vz))

    # full 3D distances between consecutive turning points
    d3_vx = np.sqrt(np.diff(px_vx)**2 + np.diff(py_vx)**2 + np.diff(pz_vx)**2)
    d3_vy = np.sqrt(np.diff(px_vy)**2 + np.diff(py_vy)**2 + np.diff(pz_vy)**2)
    d3_vz = np.sqrt(np.diff(px_vz)**2 + np.diff(py_vz)**2 + np.diff(pz_vz)**2)

    return {
        "idx_vx": idx_vx,
        "idx_vy": idx_vy,
        "idx_vz": idx_vz,

        "dx_turn": dx_turn,
        "dy_turn": dy_turn,
        "dz_turn": dz_turn,

        "d3_vx": d3_vx,
        "d3_vy": d3_vy,
        "d3_vz": d3_vz,
    }
########################################

#########
#########
# start

def plot_const_coordinate_sections(
    xs, ys, zs,
    vx, vy, vz,
    xcpxs, xcpys, xcpzs,
    dts,
    label,
    filename_prefix,
    folderforPlots=None,
    autosave=0,
    nt=None,
    num=100,
    tol=1e-3,
    figsize=(18, 12),
):
    n_col = len(dts)

    views = [
        {
            "elev": -90,
            "azim": 0,
            "filename": f"{filename_prefix}_top.png",
        },
        #{
        #    "elev": 0,
        #    "azim": 0,
        #    "filename": f"{filename_prefix}.png",
        #},
    ]

    for view in views:

        fig = plt.figure(figsize=figsize)

        for j in range(n_col):

            res = turning_point_distances(
                xs[j], ys[j], zs[j],
                vx[j], vy[j], vz[j]
            )

            #print("vertical distances from vz sign changes:")
           # print(res["dz_turn"])
            print("mean vertical distance:")
            print(np.mean(res["dz_turn"]))

            N_mod = min(len(xs[j]), len(ys[j]), len(zs[j]))
            N_cp = min(len(xcpxs[j]), len(xcpys[j]), len(xcpzs[j]))

            N = min(N_mod, N_cp)

            if nt is not None:
                N = min(N, nt)

            slice_ids = np.linspace(0, N, num, endpoint=False, dtype=int)

            modified_coords = [xs[j][:N], ys[j][:N], zs[j][:N]]
            cp_coords = [xcpxs[j][:N], xcpys[j][:N], xcpzs[j][:N]]

            row_titles = ["const. x", "const. y", "const. z"]

            for row_idx in range(3):

                ax = fig.add_subplot(
                    3,
                    n_col,
                    row_idx * n_col + 1 + j,
                    projection="3d"
                )

                ax.set_facecolor("none")

                # modified trajectory
                ax.scatter(
                    xs[j][:N],
                    ys[j][:N],
                    zs[j][:N],
                    color="magenta",
                    s=0.01,
                    marker="o"
                )

                # reference CP trajectory
                ax.plot(
                    xcpxs[j][:N],
                    xcpys[j][:N],
                    xcpzs[j][:N],
                    color="r",
                    lw=1
                )

                #ax.quiver(
                #    xs[j][::10], ys[j][::10], zs[j][::10],
                #    vx[j][::10], vy[j][::10], vz[j][::10],
                #    length=0.1,
                #    normalize=True,
                #    color="green",
                #    lw=0.1
                #)


                # coordinate of modified trajectory, e.g. x, y, or z
                coord_modified = modified_coords[row_idx]

                # coordinate of CP trajectory, e.g. x_cp, y_cp, or z_cp

                #coord_cp = modified_coords[row_idx]
                coord_cp = cp_coords[row_idx]

                for k in slice_ids:

                    # slicing value comes from CP trajectory
                    coord0 = coord_cp[k]

                    # but matching is done on the modified trajectory
                    idx = np.where(np.abs(coord_modified - coord0) < tol)[0]

                    ax.plot(
                        xs[j][idx],
                        ys[j][idx],
                        zs[j][idx],
                        color="k",
                        lw=0.1
                    )

                ax.view_init(elev=view["elev"], azim=view["azim"])

                dyj = 0.5 * dts[j]
                ax.set_ylim(-dyj, dyj)
                ax.grid(False)
                #ax.set_axis_off()
                #dxj = 10
                #ax.set_xlim(-dxj, dxj)


        col_titles = [
            rf"$\Delta t_{{\rm {label}}} = {dt}$"
            for dt in dts
        ]

        xcols = np.linspace(0.1, 0.9, n_col)

        for x, title in zip(xcols, col_titles):
            fig.text(
                x,
                0.01,
                title,
                ha="center",
                va="bottom",
                fontsize=16
            )

        row_titles = ["const. x", "const. y", "const. z"]

        for i, title in enumerate(row_titles):
            y = 1 - (i + 0.1) / 3
            fig.text(
                0.1,
                y,
                title,
                ha="right",
                va="center",
                fontsize=12
            )

        # ---- shared legend for whole figure ----
        legend_handles = [
            Line2D([0], [0], marker='o', color='none',
                   markerfacecolor='magenta', markersize=6,
                   label=f"{label} trajectory"),
            Line2D([0], [0], color='r', lw=1.5,
                   label="CP trajectory"),
            Line2D([0], [0], color='k', lw=1.5,
                   label="constant x/y/z line"),
        ]

        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.98),
            ncol=3,
            frameon=False,
            fontsize=12
        )

        plt.subplots_adjust(
            left=0.02,
            right=0.98,
            top=1,
            bottom=0,
            wspace=-0.1,
            hspace=0
        )

        if autosave == 1 and folderforPlots is not None:
            fig.savefig(
                os.path.join(folderforPlots, view["filename"]),
                dpi=300
            )

        plt.show()

#######################################################



#
# plooooooooooooooooots
#

num = 200



plot_const_coordinate_sections(
    xs=xcpx,
    ys=xcpy,
    zs=xcpz,
    vx=vcpx,
    vy=vcpy,
    vz=vcpz,
    xcpxs=xcpx,
    xcpys=xcpy,
    xcpzs=xcpz,
    dts=dt_cp_sym,
    label="cp",
    filename_prefix="const_X_cp",
    folderforPlots=folderforPlots,
    autosave=autosave,
    nt=nt,
    num=num,
    tol=1e-4,
)


plot_const_coordinate_sections(
    xs=xcppx,
    ys=xcppy,
    zs=xcppz,
    vx=vcppx,
    vy=vcppy,
    vz=vcppz,
    xcpxs=xcpx,
    xcpys=xcpy,
    xcpzs=xcpz,
    dts=dt_cpp_sym,
    label="cpp\_sym",
    filename_prefix="const_X_cppsym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    nt=nt,
    num=num,
    tol=1e-4,
)

plot_const_coordinate_sections(
    xs=xcppvarx,
    ys=xcppvary,
    zs=xcppvarz,
    vx=vcppvarx,
    vy=vcppvary,
    vz=vcppvarz,
    xcpxs=xcpx,
    xcpys=xcpy,
    xcpzs=xcpz,
    dts=dt_cpp_var,
    label="cpp\_var",
    filename_prefix="const_X_cppvar",
    folderforPlots=folderforPlots,
    autosave=autosave,
    nt=nt,
    num=num,
    tol=1e-4,
)

plot_const_coordinate_sections(
    xs=xgcx,
    ys=xgcy,
    zs=xgcz,
    vx=vgcx,
    vy=vgcy,
    vz=vgcz,
    xcpxs=xcpx,
    xcpys=xcpy,
    xcpzs=xcpz,
    dts=dt_gc_var,
    label="gc\_var",
    filename_prefix="const_X_gc",
    folderforPlots=folderforPlots,
    autosave=autosave,
    nt=nt,
    num=num,
    tol=1e-4,
)





#%%
######################## ENDE PLOTS ########################


#%%
##

from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D


def rel_plot_integrator(
    qx, qy, qz,
    vx, vy, vz,
    px, py, pz,
    dts,
    dt_label,
    integrator_name,
    filename_prefix,
    folderforPlots=None,
    autosave=0,
    q_kind="scatter",
    v_kind="plot",
    p_kind="plot",
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=0.02,
    v_weight=0.02,
    p_weight=0.02,
    figsize=(18, 18),
    other_view=None,
):
    n_col = len(dts)

    if np.isscalar(q_weight):
        q_weight = [q_weight] * n_col
    if np.isscalar(v_weight):
        v_weight = [v_weight] * n_col
    if np.isscalar(p_weight):
        p_weight = [p_weight] * n_col

    fig, axs = plt.subplots(
        3, n_col,
        figsize=figsize,
        subplot_kw={"projection": "3d"}
    )

    # in case n_col = 1, make axs always 2D-indexable
    if n_col == 1:
        axs = np.array(axs).reshape(3, 1)

    row_info = [
        {
            "name": "Position",
            "xs": qx, "ys": qy, "zs": qz,
            "kind": q_kind,
            "color": q_color,
            "weights": q_weight,
            "overlay_x": px, "overlay_y": py, "overlay_z": pz,
            "overlay_color": "r",
        },
        {
            "name": "Velocity",
            "xs": vx, "ys": vy, "zs": vz,
            "kind": v_kind,
            "color": v_color,
            "weights": v_weight,
            "overlay_x": None, "overlay_y": None, "overlay_z": None,
            "overlay_color": None,
        },
        {
            "name": "Momentum",
            "xs": px, "ys": py, "zs": pz,
            "kind": p_kind,
            "color": p_color,
            "weights": p_weight,
            "overlay_x": vx, "overlay_y": vy, "overlay_z": vz,
            "overlay_color": "g",
        },
    ]

    for i, info in enumerate(row_info):
        for j in range(n_col):
            ax = axs[i, j]

            xs = info["xs"][j]
            ys = info["ys"][j]
            zs = info["zs"][j]
            kind = info["kind"]
            color = info["color"]
            weight = info["weights"][j]


            if kind == "scatter":
                ax.scatter(xs, ys, zs, color=color, s=weight)

                if info["overlay_x"] is not None:
                    ax.scatter(
                        info["overlay_x"][j],
                        info["overlay_y"][j],
                        info["overlay_z"][j],
                        color=info["overlay_color"],
                        s=0.02
                    )

            elif kind == "plot":
                ax.plot(xs, ys, zs, color=color, lw=weight)

            else:
                raise ValueError("kind must be either 'scatter' or 'plot'")

            # topdown view
            ax.view_init(elev=-90, azim=0)
            if j == 0:
                ax.set_xlabel(r"$x$", labelpad=10, fontsize=12)
            ax.set_ylabel(r"$y$", labelpad=15, fontsize=11)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()

            ax.set_xticks(np.linspace(xmin, xmax, 4))
            ax.set_yticks(np.linspace(ymin, ymax, 4))

            #ax.xaxis.set_major_formatter(FuncFormatter(sci_fmt))
            #ax.yaxis.set_major_formatter(FuncFormatter(sci_fmt))

            ax.set_zticks([])

            ax.tick_params(axis='x', pad=20)
            ax.tick_params(axis='y', pad=10)


    handles = [
        Line2D([0], [0], color=q_color, lw=3),
        Line2D([0], [0], color=v_color, lw=3),
        Line2D([0], [0], color=p_color, lw=3),
    ]

    labels = ["Position", "Velocity", "Momentum"]

    fig.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1),   # manually adjust vertical position
        ncol=3,                       # puts legends in one row
        frameon=False,
        fontsize=15,
        handlelength=2.5,
        handletextpad=0.6,
        columnspacing=2.5
        )


    for j in range(n_col):
        y_pos = ((j + 0.5) / n_col)*0.62 + 0.2
        fig.text(
            y_pos, 0.02,
            rf"$\Delta t_{{\rm {dt_label}}} = {dts[j]}$",
            ha="center",
            va="center",
            fontsize=14
        )

    fig.subplots_adjust(
        left=0.04,
        right=0.98,
        bottom=0.04,
        top=1,
        wspace=-0.69,
        hspace=-0.102
    )

    if autosave == 1 and folderforPlots is not None:
        fig.savefig(
            os.path.join(folderforPlots, f"rel_all_{filename_prefix}.png"),
            dpi=300
        )

    plt.show()



kind = "scatter"
we = [0.05, 0.02, 0.02, 0.02]

rel_plot_integrator(
    qx=xcpx, qy=xcpy, qz=xcpz,
    vx=vcpx, vy=vcpy, vz=vcpz,
    px=pcpx, py=pcpy, pz=pcpz,
    dts=dt_cp_sym,
    dt_label="cp\_sym",
    integrator_name="Classical Particle",
    filename_prefix="cp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)




#%%

rel_plot_integrator(
    qx=xcppx, qy=xcppy, qz=xcppz,
    vx=vcppx, vy=vcppy, vz=vcppz,
    px=pcppx, py=pcppy, pz=pcppz,
    dts=dt_cpp_sym,
    dt_label="cpp\_sym",
    integrator_name="CPP Symplectic",
    filename_prefix="cpp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)

rel_plot_integrator(
    qx=xcppvarx, qy=xcppvary, qz=xcppvarz,
    vx=vcppvarx, vy=vcppvary, vz=vcppvarz,
    px=pcppvarx, py=pcppvary, pz=pcppvarz,
    dts=dt_cpp_var,
    dt_label="cpp\_var",
    integrator_name="CPP Variational",
    filename_prefix="cpp_var",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)

rel_plot_integrator(
    qx=xgcx, qy=xgcy, qz=xgcz,
    vx=vgcx, vy=vgcy, vz=vgcz,
    px=pgcx, py=pgcy, pz=pgcz,
    dts=dt_gc_var,
    dt_label="gc\_var",
    integrator_name="Guiding Center",
    filename_prefix="gc",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)



#%%
# velocity evolution, in vpar and vperp coordinates
#%%

from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D


def rel_plot_integrator(
    vparx, vpary, vparz,
    vx, vy, vz,
    vperx, vpery, vperz,
    dts,
    dt_label,
    integrator_name,
    filename_prefix,
    folderforPlots=None,
    autosave=0,
    q_kind="scatter",
    v_kind="plot",
    p_kind="plot",
    q_color="g",
    v_color="b",
    p_color="r",
    q_weight=0.02,
    v_weight=0.02,
    p_weight=0.02,
    figsize=(18, 18),
    other_view=None,
):
    n_col = len(dts)

    if np.isscalar(q_weight):
        q_weight = [q_weight] * n_col
    if np.isscalar(v_weight):
        v_weight = [v_weight] * n_col
    if np.isscalar(p_weight):
        p_weight = [p_weight] * n_col

    fig, axs = plt.subplots(
        3, n_col,
        figsize=figsize,
        subplot_kw={"projection": "3d"}
    )

    # in case n_col = 1, make axs always 2D-indexable
    if n_col == 1:
        axs = np.array(axs).reshape(3, 1)

    row_info = [
        {
            "name": "Velocity",
            "xs": vx, "ys": vy, "zs": vz,
            "kind": q_kind,
            "color": q_color,
            "weights": q_weight,
        },
        {
            "name": "Parallel Vel.",
            "xs": vparx, "ys": vpary, "zs": vparz,
            "kind": v_kind,
            "color": v_color,
            "weights": v_weight,
        },
        {
            "name": "Perpendicular Vel.",
            "xs": vperx, "ys": vpery, "zs": vperz,
            "kind": p_kind,
            "color": p_color,
            "weights": p_weight,
        },
    ]

    for i, info in enumerate(row_info):
        for j in range(n_col):
            ax = axs[i, j]

            xs = info["xs"][j]
            ys = info["ys"][j]
            zs = info["zs"][j]
            kind = info["kind"]
            color = info["color"]
            weight = info["weights"][j]


            if kind == "scatter":
                ax.scatter(xs, ys, zs, color=color, s=weight)


            elif kind == "plot":
                ax.plot(xs, ys, zs, color=color, lw=weight)

            else:
                raise ValueError("kind must be either 'scatter' or 'plot'")

            # topdown view
            ax.view_init(elev=-90, azim=0)
            if j == 0:
                ax.set_xlabel(r"$x$", labelpad=10, fontsize=12)
            ax.set_ylabel(r"$y$", labelpad=15, fontsize=11)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()

            ax.set_xticks(np.linspace(xmin, xmax, 4))
            ax.set_yticks(np.linspace(ymin, ymax, 4))

            #ax.xaxis.set_major_formatter(FuncFormatter(sci_fmt))
            #ax.yaxis.set_major_formatter(FuncFormatter(sci_fmt))

            ax.set_zticks([])

            ax.tick_params(axis='x', pad=20)
            ax.tick_params(axis='y', pad=10)


    handles = [
        Line2D([0], [0], color=q_color, lw=3),
        Line2D([0], [0], color=v_color, lw=3),
        Line2D([0], [0], color=p_color, lw=3),
    ]

    labels = ["Velocity", "Parallel Vel.", "Perpendicular Vel."]

    fig.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1),   # manually adjust vertical position
        ncol=3,                       # puts legends in one row
        frameon=False,
        fontsize=15,
        handlelength=2.5,
        handletextpad=0.6,
        columnspacing=2.5
        )


    for j in range(n_col):
        y_pos = ((j + 0.5) / n_col)*0.62 + 0.2
        fig.text(
            y_pos, 0.02,
            rf"$\Delta t_{{\rm {dt_label}}} = {dts[j]}$",
            ha="center",
            va="center",
            fontsize=14
        )

    fig.subplots_adjust(
        left=0.04,
        right=0.98,
        bottom=0.04,
        top=1,
        wspace=-0.69,
        hspace=-0.102
    )

    if autosave == 1 and folderforPlots is not None:
        fig.savefig(
            os.path.join(folderforPlots, f"rel_allVEL_{filename_prefix}.png"),
            dpi=300
        )
    plt.show()




kind = "scatter"
we = [0.05, 0.02, 0.02, 0.02]

rel_plot_integrator(
    vparx=vparcpx, vpary=vparcpy, vparz=vparcpz,
    vx=vcpx, vy=vcpy, vz=vcpz,
    vperx=vpercpx, vpery=vpercpy, vperz=vpercpz,
    dts=dt_cpp_sym,
    dt_label="cp\_sym",
    integrator_name="CP Symplectic",
    filename_prefix="cp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="seagreen",
    v_color="steelblue",
    p_color="darkorange",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)


#%%

rel_plot_integrator(
    vparx=vparcppx, vpary=vparcppy, vparz=vparcppz,
    vx=vcppx, vy=vcppy, vz=vcppz,
    vperx=vpercppx, vpery=vpercppy, vperz=vpercppz,
    dts=dt_cpp_sym,
    dt_label="cpp\_sym",
    integrator_name="CPP Symplectic",
    filename_prefix="cpp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="seagreen",
    v_color="steelblue",
    p_color="darkorange",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)

rel_plot_integrator(
    vparx=vparcppvarx, vpary=vparcppvary, vparz=vparcppvarz,
    vx=vcppvarx, vy=vcppvary, vz=vcppvarz,
    vperx=vpercppvarx, vpery=vpercppvary, vperz=vpercppvarz,
    dts=dt_cpp_var,
    dt_label="cpp\_var",
    integrator_name="CPP Variational",
    filename_prefix="cpp_var",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="seagreen",
    v_color="steelblue",
    p_color="darkorange",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)

rel_plot_integrator(
    vparx=vpargcx, vpary=vpargcy, vparz=vpargcz,
    vx=vgcx, vy=vgcy, vz=vgcz,
    vperx=vpergcx, vpery=vpergcy, vperz=vpergcz,
    dts=dt_gc_var,
    dt_label="gc\_var",
    integrator_name="Guiding Center",
    filename_prefix="gc",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="seagreen",
    v_color="steelblue",
    p_color="darkorange",
    q_weight=we,
    v_weight=we,
    p_weight=we,
    figsize=(22, 12)
)



#%%

#
#
# ###
# Position --> Blue
# ###
# Velocity --> Green
# ###
# Momentum --> Red
# ###
#
#
#

from matplotlib.ticker import FuncFormatter

def sci_fmt(x, pos):
    if np.isclose(x, 0):
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    coeff = np.round(x / 10**exp, 3)
    return rf"${np.round(coeff):g}\cdot 10^{{{exp}}}$"



def rel_plot_integrator(
    qx, qy, qz,
    vx, vy, vz,
    px, py, pz,
    dts,
    dt_label,
    integrator_name,
    filename_prefix,
    folderforPlots=None,
    autosave=0,
    q_kind="scatter",
    v_kind="plot",
    p_kind="plot",
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=0.02,
    v_weight=0.02,
    p_weight=0.02,
    figsize_q=(18, 6),
    figsize_v=(18, 6),
    figsize_p=(18, 6),
    other_view_q=None,
    other_view_v=None,
    other_view_p=None,
):
    """
    For one integrator, this creates four figures:

        1. position topdown
        2. position other view
        3. velocity topdown
        4. velocity other view
        5. momentum topdown
        6. momentum other view

    Example filename_prefix:
        'cp'
        'cpp_sym'
        'cpp_var'
        'gc'
    """

    n_col = len(dts)

    if np.isscalar(q_weight):
        q_weight = [q_weight] * n_col

    if np.isscalar(v_weight):
        v_weight = [v_weight] * n_col

    if np.isscalar(p_weight):
        p_weight = [p_weight] * n_col

    views = [
        {
            "name": "top",
            "suffix": "_top",
            "title_suffix": " - Topdown",
            "elev": -90,
            "azim": 0,
        },
        #{
        #    "name": "other",
        #    "suffix": "",
        #    "title_suffix": " - other view",
        #    "elev": None,
        #    "azim": None,
        #},
    ]

    for quantity in ["position", "velocity", "momentum"]:

        if quantity == "position":
            xs, ys, zs = qx, qy, qz
            kind = q_kind
            color = q_color
            weights = q_weight
            figsize = figsize_q
            suptitle_base = f"Relative position {integrator_name}"
            filename_base = f"rel_q{filename_prefix}"
            other_view = other_view_q

        elif quantity == "velocity":
            xs, ys, zs = vx, vy, vz
            kind = v_kind
            color = v_color
            weights = v_weight
            figsize = figsize_v
            suptitle_base = f"Velocity of {integrator_name}"
            filename_base = f"rel_v{filename_prefix}"
            other_view = other_view_v

        elif quantity == "momentum":
            xs, ys, zs = px, py, pz
            kind = p_kind
            color = p_color
            weights = p_weight
            figsize = figsize_p
            suptitle_base = f"Momentum of {integrator_name}"
            filename_base = f"rel_p{filename_prefix}"
            other_view = other_view_p

        for view in views:

            fig = plt.figure(figsize=figsize)

            for j in range(n_col):

                ax = fig.add_subplot(
                    1,
                    n_col,
                    j + 1,
                    projection="3d"
                )

                if kind == "scatter":
                    ax.scatter(
                        xs[j],
                        ys[j],
                        zs[j],
                        color=color,
                        s=weights[j]
                    )
                    if quantity == "position":

                        ax.scatter(
                            px[j],
                            py[j],
                            pz[j],
                            color = 'r',
                            s = 0.02)

                    elif quantity == "momentum":

                        ax.scatter(
                            vx[j],
                            vy[j],
                            vz[j],
                            color = 'g',
                            s = 0.02)

                elif kind == "plot":
                    ax.plot(
                        xs[j],
                        ys[j],
                        zs[j],
                        color=color,
                        lw=weights[j]
                    )

                else:
                    raise ValueError("kind must be either 'scatter' or 'plot'")

                ax.set_title(
                    rf"$\Delta t_{{\rm {dt_label}}} = {dts[j]}$",
                    fontsize=14
                )

                if view["name"] == "top":
                    ax.view_init(elev=view["elev"], azim=view["azim"])

                elif view["name"] == "other" and other_view is not None:
                    ax.view_init(elev=other_view[0], azim=other_view[1])



                ax.set_xlabel(r"$x$", labelpad=10, fontsize=11)
                ax.set_ylabel(r"$y$", labelpad=10, fontsize=11)


                xmin, xmax = ax.get_xlim()
                ymin, ymax = ax.get_ylim()

                ax.set_xticks(np.linspace(xmin, xmax, 4))
                ax.set_yticks(np.linspace(ymin, ymax, 4))
                ax.xaxis.set_major_formatter(FuncFormatter(sci_fmt))
                ax.yaxis.set_major_formatter(FuncFormatter(sci_fmt))
                ax.set_zticks([])

                ax.tick_params(axis='x', pad=20)

            fig.subplots_adjust(left=0.03, right=0.98, bottom=0., top=1, wspace=-0.2)

            if autosave == 1 and folderforPlots is not None:

                fig.savefig(
                    os.path.join(
                        folderforPlots,
                        filename_base + view["suffix"] + ".png"
                    ),
                    dpi=300
                )

            plt.show()




kind = "scatter"
we = [0.05, 0.02, 0.02, 0.02]

# scale: vel -> q : 1e5 * vel = pos
#        vel -> p : 30 * vel =  ((2*pi)^2 ~ 39 ?)
#        p -> q : 1e5 / 30 ~ 3333
# Classical Particle Symplectic

rel_plot_integrator(
    qx=xcpx, qy=xcpy, qz=xcpz,
    vx=vcpx, vy=vcpy, vz=vcpz,
    px=pcpx, py=pcpy, pz=pcpz,
    dts=dt_cp_sym,
    dt_label="cp\_sym",
    integrator_name="Classical Particle",
    filename_prefix="cp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we
)


#%%

# Classical Pauli Particle Symplectic

rel_plot_integrator(
    qx=xcppx, qy=xcppy, qz=xcppz,
    vx=vcppx, vy=vcppy, vz=vcppz,
    px=pcppx, py=pcppy, pz=pcppz,
    dts=dt_cpp_sym,
    dt_label="cpp\_sym",
    integrator_name="CPP symplectic",
    filename_prefix="cpp_sym",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we
)


#%%

# Classical Pauli Particle Variational

rel_plot_integrator(
    qx=xcppvarx, qy=xcppvary, qz=xcppvarz,
    vx=vcppvarx, vy=vcppvary, vz=vcppvarz,
    px=pcppvarx, py=pcppvary, pz=pcppvarz,
    dts=dt_cpp_var,
    dt_label="cpp\_var",
    integrator_name="CPP variational",
    filename_prefix="cpp_var",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we
)


# %%

rel_plot_integrator(
    qx=xgcx, qy=xgcy, qz=xgcz,
    vx=vgcx, vy=vgcy, vz=vgcz,
    px=pgcx, py=pgcy, pz=pgcz,
    dts=dt_gc_var,
    dt_label="gc,var",
    integrator_name="Guiding Center",
    filename_prefix="gc_var",
    folderforPlots=folderforPlots,
    autosave=autosave,
    q_kind=kind,
    v_kind=kind,
    p_kind=kind,
    q_color="b",
    v_color="green",
    p_color="r",
    q_weight=we,
    v_weight=we,
    p_weight=we
)

# %%



def plot_xyz_histograms(
    xs, ys, zs,
    dts,
    label,
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=None,
    bins=100,
    figsize=None,
    color=None
):
    """
    Plot histograms of x, y, z coordinates.

    Rows:
        1. x-coordinate histogram
        2. y-coordinate histogram
        3. z-coordinate histogram

    Columns:
        different time steps in dts
    """

    n_col = len(dts)

    if figsize is None:
        figsize = (4.5 * n_col, 9)

    fig, axes = plt.subplots(3, n_col, figsize=figsize)

    # Make axes always two-dimensional, also for n_col = 1
    if n_col == 1:
        axes = axes.reshape(3, 1)

    data_list = [xs, ys, zs]
    coord_labels = [r"$x$", r"$y$", r"$z$"]

    for i in range(3):
        data = data_list[i]

        for j in range(n_col):

            N = len(data[j])

            if nt is not None:
                N = min(N, nt)

            ax = axes[i, j]

            ax.hist(
                data[j][:N],
                bins=bins,
                color=color,
                edgecolor="k"
            )

            ax.grid(True)

            # Column title only in the first row
            if i == 0:
                ax.set_title(rf"$\Delta t_{{\rm {label}}} = {dts[j]}$")

            ax.set_xlabel(coord_labels[i])
            ax.set_ylabel("counts")

    plt.tight_layout()

    if autosave == 1 and folderforPlots is not None:
        fig.savefig(
            os.path.join(folderforPlots, f"hist_xyz_{label}.png"),
            dpi=300,
            bbox_inches="tight"
        )

    plt.show()


plot_xyz_histograms(
    xs=xcpx,
    ys=xcpy,
    zs=xcpz,
    dts=dt_cp_sym,
    label="xcp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightblue"
)
plot_xyz_histograms(
    xs=vcpx,
    ys=vcpy,
    zs=vcpz,
    dts=dt_cp_sym,
    label="vcp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightgreen"
)
plot_xyz_histograms(
    xs=pcpx,
    ys=pcpy,
    zs=pcpz,
    dts=dt_cp_sym,
    label="pcp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="salmon"
)




#%%
# CPP Symplectic

plot_xyz_histograms(
    xs=xcppx,
    ys=xcppy,
    zs=xcppz,
    dts=dt_cpp_sym,
    label="xcpp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightblue"
)
plot_xyz_histograms(
    xs=vcppx,
    ys=vcppy,
    zs=vcppz,
    dts=dt_cpp_sym,
    label="vcpp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightgreen"
)
plot_xyz_histograms(
    xs=pcppx,
    ys=pcppy,
    zs=pcppz,
    dts=dt_cpp_sym,
    label="pcpp\_sym",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="salmon"
)

#%%
# CPP Variational

plot_xyz_histograms(
    xs=xcppvarx,
    ys=xcppvary,
    zs=xcppvarz,
    dts=dt_cpp_var,
    label="xcpp\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightblue"
)
plot_xyz_histograms(
    xs=vcppvarx,
    ys=vcppvary,
    zs=vcppvarz,
    dts=dt_cpp_var,
    label="vcpp\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightgreen"
)
plot_xyz_histograms(
    xs=pcppvarx,
    ys=pcppvary,
    zs=pcppvarz,
    dts=dt_cpp_var,
    label="pcpp\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="salmon"
)

#%%

plot_xyz_histograms(
    xs=xgcx,
    ys=xgcy,
    zs=xgcz,
    dts=dt_gc_var,
    label="xgc\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightblue"
)
plot_xyz_histograms(
    xs=vgcx,
    ys=vgcy,
    zs=vgcz,
    dts=dt_gc_var,
    label="vgc\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="lightgreen"
)
plot_xyz_histograms(
    xs=pgcx,
    ys=pgcy,
    zs=pgcz,
    dts=dt_gc_var,
    label="pgc\_var",
    autosave=autosave,
    folderforPlots=folderforPlots,
    nt=nt,
    bins=100,
    color="salmon"
)

# %%
