# %%
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
from plotale import plot_orbit, plot_mani, plot_cost_function
from cpp_sym_midpoint import z, vel
from cp_sym_midpoint import q_cp, v_cp

'''
def position(M, z0):

    distances = np.linalg.norm(M - z0, axis=1)

    # Find index of the minimum distance
    idx_min = np.argmin(distances)
'''
traj_v = vel
traj_q = z

traj_v_cp = v_cp
traj_q_cp = q_cp    

f = field()

R0 = 1
z0 = 1
print(r0, th0, ph0)

#Resolution
n = 100#200 # not to be odd !!
nh = 100#100 # n/2
l = n**2 * nh # Number of Pixels



# Grid in flux coordinates

gridR = np.linspace(0,0.5,n)
gridTh = np.linspace(0,2*np.pi,n)
gridPh = np.linspace(0,2*np.pi,n)

a = 0.5
th, ph = np.meshgrid(gridTh, gridPh, indexing='ij')

X = (R0 + a*np.cos(th))*np.cos(ph)
Y = (R0 + a*np.cos(th))*np.sin(ph)
Z = z0 + a*np.sin(th)



#fig = plt.figure()

#ax = fig.add_subplot(111, projection='3d')

#sc = ax.scatter(X, Y, Z, cmap='viridis', s=0.1)

#plt.show()




# Cartesian grid
# Components x^i
X = np.linspace(-3,3,n)
Y = np.linspace(-3,3,n)
Z = np.linspace(0,2,nh)

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
                #elif x > 0 and y < 0:
                #    ph = 2*np.pi + phi
                else: break

                q[c,:] = [r, th, ph]
                c +=1

q = q[:c,:]
print(c)

# x^i = x^i(q)
x = (R0 + q[:,0] * np.cos(q[:,1])) * np.cos(q[:,2])
y = (R0 + q[:,0] * np.cos(q[:,1])) * np.sin(q[:,2])
z = z0 + q[:,0] * np.sin(q[:,1])
###
#
#
#
#
#
#
#
#


    #vx5[j] = x_dot[0]
    #vy5[j] = x_dot[1]
    #vz5[j] = x_dot[2]



#velx5 = (R0 + traj_v[0,:500] * np.cos(traj_v[1,:500])) * np.cos(traj_v[2,:500])
#vely5 = (R0 + traj_v[0,:500] * np.cos(traj_v[1,:500])) * np.sin(traj_v[2,:500])
#velz5 = z0 + traj_v[0,:500] * np.sin(traj_v[1,:500]) 

# cp:

vcp_x = np.zeros([3,len(traj_q_cp[0,:])])

for j in range(len(traj_q_cp[0,:])):

    A = np.array([ [np.cos(traj_q_cp[1,j])*np.cos(traj_q_cp[2,j]), -r*np.sin(traj_q_cp[1,j])*np.cos(traj_q_cp[2,j]), -np.sin(traj_q_cp[2,j])], 
                      [np.cos(traj_q_cp[1,j])*np.sin(traj_q_cp[2,j]), -r*np.sin(traj_q_cp[1,j])*np.sin(traj_q_cp[2,j]),  np.cos(traj_q_cp[2,j])],
                      [               np.sin(traj_q_cp[1,j]),                 r*np.cos(traj_q_cp[1,j]),                 0] ])

    q_dot = traj_v_cp[:,j]

    x_dot = A @ q_dot

    vcp_x[:,j] = x_dot


# Trajectory of particle in Cartesian coordinates

qcpx = (R0 + traj_q_cp[0,:] * np.cos(traj_q_cp[1,:])) * np.cos(traj_q_cp[2,:])
qcpy = (R0 + traj_q_cp[0,:] * np.cos(traj_q_cp[1,:])) * np.sin(traj_q_cp[2,:])
qcpz = z0 + traj_q_cp[0,:] * np.sin(traj_q_cp[1,:])


# Velocity in Cartesian coordinates
vcpx = (vcp_x[0,:]) # + qcpz
vcpy = (vcp_x[1,:]) # + qcpy
vcpz = (vcp_x[2,:]) # + qcpz



# cpp:

# v^x = dx/dt = dx/dq * vel^q
vel_x = np.zeros([3,len(traj_q[0,:])])

for j in range(len(traj_q[0,:])):

    A = np.array([ [np.cos(traj_q[1,j])*np.cos(traj_q[2,j]), -r*np.sin(traj_q[1,j])*np.cos(traj_q[2,j]), -np.sin(traj_q[2,j])], 
                      [np.cos(traj_q[1,j])*np.sin(traj_q[2,j]), -r*np.sin(traj_q[1,j])*np.sin(traj_q[2,j]),  np.cos(traj_q[2,j])],
                      [               np.sin(traj_q[1,j]),                 r*np.cos(traj_q[1,j]),                 0] ])

    q_dot = traj_v[:,j]

    x_dot = A @ q_dot

    vel_x[:,j] = x_dot

# Trajectory of particle in Cartesian coordinates

z1 = (R0 + traj_q[0,:] * np.cos(traj_q[1,:])) * np.cos(traj_q[2,:])
z2 = (R0 + traj_q[0,:] * np.cos(traj_q[1,:])) * np.sin(traj_q[2,:])
z3 = z0 + traj_q[0,:] * np.sin(traj_q[1,:])


# first 50 evolutions of the trajectory

x50 = (R0 + traj_q[0,:50] * np.cos(traj_q[1,:50])) * np.cos(traj_q[2,:50])
y50 = (R0 + traj_q[0,:50] * np.cos(traj_q[1,:50])) * np.sin(traj_q[2,:50])
z50 = z0 + traj_q[0,:50] * np.sin(traj_q[1,:50])


# Velocity in Cartesian coordinates
velx = (vel_x[0,:])#*50  + z1
vely = (vel_x[1,:])#*50  + z2
velz = (vel_x[2,:])#*50  + z3






# Calculate Surfaces and Curves at Point p0
#
#
#
r0, th0, ph0 = traj_q[0,:50], traj_q[1,:50], traj_q[2,:50]
#
#
#
surf_R = []
surf_Th = []
surf_Phi = []
curv_R = []
curv_Th = []
curv_Phi = []

for i in range(0,50):
    print(i)
    r0i = r0[i]
    th0i = th0[i]
    ph0i = ph0[i]

    surf_Ri = np.c_[ (R0 + r0i*np.cos(q[::100,1]))*np.cos(q[::100,2]),
                     (R0 + r0i*np.cos(q[::100,1]))*np.sin(q[::100,2]),
                     z0 + r0i*np.sin(q[::100,1]) ]

    surf_Thi = np.c_[ (R0 + q[:,0]*np.cos(th0i))*np.cos(q[:,2]),
                      (R0 + q[:,0]*np.cos(th0i))*np.sin(q[:,2]),
                      z0 + q[:,0]*np.sin(th0i) ]

    surf_Phii = np.c_[ (R0 + q[:,0]*np.cos(q[:,1]))*np.cos(ph0i),
                       (R0 + q[:,0]*np.cos(q[:,1]))*np.sin(ph0i),
                       z0 + q[:,0]*np.sin(q[:,1]) ]

    curv_Ri = np.c_[ (R0 + q[::100,0]*np.cos(th0i))*np.cos(ph0i),
                      (R0 + q[::100,0]*np.cos(th0i))*np.sin(ph0i),
                      z0 + q[::100,0]*np.sin(th0i) ]

    curv_Thi = np.c_[ (R0 + r0i*np.cos(q[::100,1]))*np.cos(ph0i),
                       (R0 + r0i*np.cos(q[::100,1]))*np.sin(ph0i),
                       z0 + r0i*np.sin(q[::100,1]) ]

    curv_Phii = np.c_[ (R0 + r0i*np.cos(th0i))*np.cos(q[::100,2]),
                       (R0 + r0i*np.cos(th0i))*np.sin(q[::100,2]),
                       z0 + r0i*np.sin(th0i)*np.ones(np.size(q[::100,2])) ]

    surf_R.append(surf_Ri)
    surf_Th.append(surf_Thi)
    surf_Phi.append(surf_Phii)
    curv_R.append(curv_Ri)
    curv_Th.append(curv_Thi)
    curv_Phi.append(curv_Phii)



# Interpolation of the field on the grid points


# Field evaluation on trajectory at grid points 

B = np.zeros(c)
k = 0
zero = 1000


for i in range(0,c):
    r = q[i,0]
    th = q[i,1]
    ph = q[i,2]
    f.evaluate(r,th,ph)
    B[i] = f.B


    if i == 0:
        const_B = []
        f.evaluate(q[zero,0], q[zero,1], q[zero,2])
        const_B.append([q[zero,0], q[zero,1], q[zero,2], f.B])
        eps = abs(const_B[i][3] - B[0])
        print(eps)
    
    elif eps <= 1e-4:
        const_B.append([q[i-1,0], q[i-1,1], q[i-1,2], B[i-1]])
        eps = abs(B[i] - const_B[0][3]) 
        k += 1

    else:
        eps = abs(const_B[0][3] - B[i])
    




print('c =' , c)
print('k =' , k)
# with k the number of points with almost the same B field value
# and c the total number of grid points inside the torus
arr = np.array(const_B)
print(np.shape(arr))
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


# plot Banana orbit (Curie plot idk)

fig, ax = plt.subplots()

plot_orbit(traj_q, ax=ax)
plot_orbit(traj_q_cp, ax=ax)
plot_orbit(consB, ax=ax)
plt.show()


fig, ax = plt.subplots()

plot_orbit(traj_q, ax=ax)
plot_orbit(traj_q_cp, ax=ax)
plt.show()



# cp: Plot evverything

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(qcpx, qcpy, qcpz, c='r', cmap='viridis', s=1)
sc = ax.scatter(vcpx+qcpx, vcpy+qcpy, vcpz+qcpz, color='green', marker='o', s=1, label="Velocity points")
#sc = ax.scatter(vcpx, vcpy, vcpz, color='green', marker='o', s=1, label="Velocity relative")
for (q1, th1, ph1, q2, th2, ph2) in zip(qcpx, qcpy, qcpz, qcpx+vcpx, qcpy+vcpy, qcpz+vcpz):
    ax.plot([q1, q2], [th1, th2], [ph1, ph2],
            color='gray', linestyle='--', linewidth=0.8)

sc = ax.plot(qcpx[:50]+vcpx[:50], qcpy[:50]+vcpy[:50], qcpz[:50]+vcpz[:50], color='g', linewidth=1.5, label="Line through 50 points")

plt.title('Classical Particle in Configuration Space Q')
#ax.set_xlim([0.95,1.05])
#ax.set_ylim([-0.1,0])
#ax.set_zlim([1.005,1.105])
plt.show()




# cpp: Plot everything

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')


#sc = ax.scatter(x, y, z, c=B, cmap='viridis', s=0.1)

#sc = ax.scatter(consBx, consBy, consBz, c=B_const, cmap='viridis', marker='o', s=50, label="Constant B")
sc = ax.plot(consBx, consBy, consBz, c='b', linewidth=1.5, label="Constant B")

#sc = ax.scatter(surf_R[10][:,0], surf_R[10][:,1], surf_R[10][:,2], color='b', marker='o', s=20)#marker='o', s=200, cmap='viridis')
#sc = ax.scatter(surf_Th[10][:,0], surf_Th[10][:,1], surf_Th[10][:,2], c=B, cmap='viridis', s=0.1)
#sc = ax.scatter(surf_Phi[:,0], surf_Phi[:,1], surf_Phi[:,2], c=B, cmap='viridis', s=0.1)

#sc = ax.scatter(curv_R[10][:,0], curv_R[10][:,1], curv_R[10][:,2], color='g', s=0.5)
#sc = ax.scatter(curv_Th[10][:,0], curv_Th[10][:,1], curv_Th[10][:,2], color='g', s=0.5)
#sc = ax.scatter(curv_Phi[10][:,0], curv_Phi[10][:,1], curv_Phi[10][:,2], color='g', s=0.5)



### Trajectory of the particle: zi(t) = xi(q(t)) = f(T(t)), with q(t) being the curves on config manif Q. 
# eqiv: zi(t) \in R <=> R -> Q -> R  ---> 'Every value of 't' maps with a curve to a unique point in Q' ---> the coord. func. xi(T(t)) 

sc = ax.scatter(z1, z2, z3, c='r', cmap='viridis', s=1)

#ax.plot(x50, y50, z50, color='red', linewidth=1.5, label="Line through 50 points")
#ax.scatter(x50[::10], y50[::10], z50[::10], color='red', marker='o')
#sc = ax.scatter(x50, y50, z50, c='red', cmap='viridis', s=1)

### Plot velocity at each point of the trajectory

# Scatter of velocity points
sc = ax.scatter(velx+z1, vely+z2, velz+z3, color='green', marker='o', s=1, label="Velocity points")
#sc = ax.scatter(velx, vely, velz, color='green', marker='o', s=1, label="Velocity relative")

# Connect corresponding points from (x,y,z) to (velx,vely,velz)
for (x1, x2, x3, v1, v2, v3) in zip(z1, z2, z3, velx*50+z1, vely*50+z2, velz*50+z3):
    ax.plot([x1, v1], [x2, v2], [x3, v3],
            color='gray', linestyle='--', linewidth=0.8)


for (x1, y1, z1, x2, y2, z2) in zip(z1[::50], z2[::50], z3[::50], velx[::50], vely[::50], velz[::50]):
    ax.plot([x1, x2], [y1, y2], [z1, z2],
            color='gray', linestyle='--', linewidth=0.8)



#ax.set_xlim([-1,1])
#ax.set_ylim([0,2])
#ax.set_zlim([0,2])
plt.show()


