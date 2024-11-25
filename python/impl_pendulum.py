import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from numpy import array, zeros, arange
from scipy.optimize import root
from scipy.interpolate import lagrange
import common
from common import f, r0, th0, ph0, pph0, timesteps, get_val, get_der
from plotting import plot_orbit, plot_cost_function






L = 1
m = 1
g = 9.81
dt = 1e-2
phy0 = -np.pi/2
pphy0 = 0

N = 1000
time = np.linspace(0, dt*N, N)

#explicit
q = np.zeros(N)
q[0] = phy0
p = np.zeros(N)
p[0] = pphy0

for t in range(N-1):
    q[t+1] = q[t] + dt*(p[t] / (m*L**2))
    p[t+1] = p[t] - dt*(m*g*L*np.sin(q[t]))



#implicit
def newt(wn, mn):
    q = wn

    f = lambda x: x + dt**2 *g *np.sin(x) - dt/(m*L**2) *mn - wn
    df = lambda x: dt**2 *g *np.cos(x) + 1
    e = abs(f(q))
 
    while e >= 1e-4: 
        qnew = q - f(q)/df(q)
        # qnew = q - np.linalg.solve(J, f(q))
        q = qnew
        e = abs(f(q))    
    return q


w = np.zeros(N)
w[0] = phy0
v = np.zeros(N)
v[0] = pphy0

for t in range(N-1):
    qnew = newt(w[t], v[t])
    w[t+1] = qnew
    v[t+1] = v[t] - dt *m *g *L *np.sin(qnew)




x = newt(q[39],p[39])

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 8))  # 4 rows, 1 column

# First subplot
ax1.plot(time, q)
ax1.set_title('position')



# Second subplot
ax2.plot(time, p)
ax2.set_title('momentum')

# Third subplot
ax3.plot(time, w)
ax3.set_title('implicit_q')

# 4th subplot
ax4.plot(time, v)
ax4.set_title('implicit_p')

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()

