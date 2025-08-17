# %%
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import (r0,th0,ph0,pph0,timesteps,get_val,get_der,mu,qe,m,c,)
from field_correct_test import field
#from cpp_sym_midpoint import F

'''
def position(M, z0):

    distances = np.linalg.norm(M - z0, axis=1)

    # Find index of the minimum distance
    idx_min = np.argmin(distances)
'''


f = field()

R0 = 1
z0 = 1

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

X = np.linspace(-2,2,n)
Y = np.linspace(-2,2,n)
Z = np.linspace(0,2,nh)

# Transformation into flux coordinates 
P = np.zeros([l,3])
c = 0

for x in X:
    for y in Y:
        count = -1
        R = np.sqrt((R0 - np.sqrt(x**2 + y**2))**2 + (Z - z0)**2)

        for z in Z:
            count += 1
            p = np.zeros(3)
            
            if R[count] <= 0.5:
                
                p[0] = R[count]

                a = z-z0
                b = np.sqrt(p[0]**2 - a**2)
                theta = np.arctan(a/b)
                cond = np.sqrt(x**2 + y**2) - R0
    
                if a >= 0 and cond >= 0:
                    p[1] = theta
                elif a >= 0 and cond < 0: 
                    p[1] = np.pi - theta
                elif a < 0 and cond < 0:
                    p[1] = np.pi - theta
                elif a < 0 and cond >= 0:
                    p[1] = 2*np.pi + theta
           

                phi = np.arctan(y/x)
                if x > 0 and y > 0:
                    p[2] = phi
                elif x < 0 and y > 0: 
                    p[2] = np.pi + phi
                #elif x < 0 and y < 0:
                #    p[2] = np.pi + phi
                #elif x > 0 and y < 0:
                #    p[2] = 2*np.pi + phi
                else: break

                P[c,:] = p
                c +=1


x = (R0 + P[:c,0] * np.cos(P[:c,1])) * np.cos(P[:c,2])
y = (R0 + P[:c,0] * np.cos(P[:c,1])) * np.sin(P[:c,2])
z = z0 + P[:c,0] * np.sin(P[:c,1])

# Calculate Surfaces and Curves at Point k 
k = 20000
#print(x[k])
#print(y[k])
#print(z[k])


sur_rx = (R0 + P[k,0] * np.cos(P[:c,1])) * np.cos(P[:c,2])
sur_ry = (R0 + P[k,0] * np.cos(P[:c,1])) * np.sin(P[:c,2])
sur_rz = z0 + P[k,0] * np.sin(P[:c,1])

surf_R = np.column_stack((sur_rx, sur_ry, sur_rz))


sur_thx = (R0 + P[:c,0] * np.cos(P[k,1])) * np.cos(P[:c,2])
sur_thy = (R0 + P[:c,0] * np.cos(P[k,1])) * np.sin(P[:c,2])
sur_thz = z0 + P[:c,0] * np.sin(P[k,1])
# since r*cos(pi) -> -r*cos(0) but i also want to plot surface over full cross section
neg_sur_thx = (R0 - P[:c,0] * np.cos(P[k,1])) * np.cos(P[:c,2])
neg_sur_thy = (R0 - P[:c,0] * np.cos(P[k,1])) * np.sin(P[:c,2])
neg_sur_thz = z0 - P[:c,0] * np.sin(P[k,1])

sur_phx = (R0 + P[:c,0] * np.cos(P[:c,1])) * np.cos(P[k,2])
sur_phy = (R0 + P[:c,0] * np.cos(P[:c,1])) * np.sin(P[k,2])
sur_phz = z0 + P[:c,0] * np.sin(P[:c,1])


cur_rx = (R0 + P[:c,0] * np.cos(P[k,1])) * np.cos(P[k,2])
cur_ry = (R0 + P[:c,0] * np.cos(P[k,1])) * np.sin(P[k,2])
cur_rz = z0 + P[:c,0] * np.sin(P[k,1])

cur_thx = (R0 + P[k,0] * np.cos(P[:c,1])) * np.cos(P[k,2])
cur_thy = (R0 + P[k,0] * np.cos(P[:c,1])) * np.sin(P[k,2])
cur_thz = z0 + P[k,0] * np.sin(P[:c,1])

cur_phx = (R0 + P[k,0] * np.cos(P[k,1])) * np.cos(P[:c,2])
cur_phy = (R0 + P[k,0] * np.cos(P[k,1])) * np.sin(P[:c,2])
cur_phz = z0 + P[k,0] * np.sin(P[k,1])

B = np.zeros(c)

for i in range(0,c):
    
    r = P[i,0]
    th = P[i,1]
    ph = P[i,2]
    f.evaluate(r,th,ph)
    B[i] = f.B


#fig = plt.figure()

#ax = fig.add_subplot(111, projection='3d')

#sc = ax.scatter(x, y, z, c=B, cmap='viridis', s=0.1)
#sc = ax.scatter(sur_rx, sur_ry, sur_rz, c=B, cmap='viridis', s=0.1)
#sc = ax.scatter(sur_thx, sur_thy, sur_thz, c=B, cmap='viridis', s=0.1)
#sc = ax.scatter(neg_sur_thx, neg_sur_thy, neg_sur_thz, c=B, cmap='viridis', s=0.1)
#sc = ax.scatter(sur_phx, sur_phy, sur_phz, c=B, cmap='viridis', s=0.1)
#ax.scatter(x[k], y[k], z[k], color='red' )
#sc = ax.scatter(cur_rx, cur_ry, cur_rz, color='r', s=0.5)
#sc = ax.scatter(cur_thx, cur_thy, cur_thz, color='r', s=0.5)
#sc = ax.scatter(cur_phx, cur_phy, cur_phz, color='r', s=0.5)


#plt.show()




'''


#Surface on which phi =! 0:

indices = [i for i, val in enumerate(P[:c,2]) if val < 2e-2]
print(indices)
print(np.size(indices))

proX = np.zeros(np.size(indices))
proY = np.zeros(np.size(indices))
proZ = np.zeros(np.size(indices))
proB = np.zeros(np.size(indices))
c = 0
for i in indices:
    proX[c] = (R0 + P[i,0] * np.cos(P[i,1])) * np.cos(P[i,2])
    proY[c] = (R0 + P[i,0] * np.cos(P[i,1])) * np.sin(P[i,2])
    proZ[c] = z0 + P[i,0] * np.sin(P[i,1])
    f.evaluate(P[i,0], P[i,1], P[i,2])
    B[c] = f.B
    c += 1

#plotting basis
x1 = proX[1400]
y1 = proY[1400]
z1 = proZ[1400]

r1 = P[indices[1400],0]
th1 = P[indices[1400],1]
ph1 = P[indices[1400],2]

r2 = r1 + 0.5
th2 = np.linspace(th1, th1+np.pi, 20)
#th2 = th1 + np.pi/4
ph2 = np.linspace(ph1, ph1+np.pi*2, 20)
#ph2 = ph1 + np.pi/4

x2r = (R0 + r2*np.cos(th1)) * np.cos(ph1)
y2r = (R0 + r2*np.cos(th1)) * np.sin(ph1)
z2r = z0 + r2 *np.sin(th1)

x2th = (R0 + r1*np.cos(th2)) * np.cos(ph1)
y2th = (R0 + r1*np.cos(th2)) * np.sin(ph1)
z2th = z0 + r1 *np.sin(th2)

x2ph = (R0 + r1*np.cos(th1)) * np.cos(ph2)
y2ph = (R0 + r1*np.cos(th1)) * np.sin(ph2)
z2ph = z0 + r1 *np.sin(th1)


#x_r = [x1, x2r]
#y_r = [y1, y2r]
#z_r = [z1, z2r]
x_r = [R0, x1, x2r]
y_r = [0, y1, y2r]
z_r = [z0, z1, z2r]
x_th = x2th

y_th = y2th
z_th = z2th
x_ph = x2ph
y_ph = y2ph
z_ph = z2ph

print(x_th)

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(proX, proY, proZ, c=proB, cmap='viridis', s=0.1)
ax.scatter(proX[1400], proY[1400], proZ[1400], color='red')

ax.plot(x_r, y_r, z_r, color='g', linewidth=1)
ax.plot(x_th, y_th, z_th, color='red', linewidth=1)
ax.plot(x_ph, y_ph, z_ph, color='b', linewidth=1)



'''