""" 
Coupled equations.
"""

from dolfin import *
import numpy as np
from RKSolver import *
from pc_solver import *
import pylab as plt
        
mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, 3, 3)
V = FunctionSpace(mesh, "CG", 1)

# Unknowns
y1 = Function(V)
y2 = Function(V)

# Initial conditions
y1.vector()[:] = 3.0
y2.vector()[:] = 5.0

# Define LV parameters
A  = 1.0
B  = 1.0
C  = 1.0
D  = 1.0

# Slope functions
def fy1(t, ys):
  y1 = ys[0]
  y2 = ys[1]
  return (A * y1 - B * y1 * y2)

def fy2(t, ys):
  y1 = ys[0]
  y2 = ys[1]
  return (-C * y2 + D * y1 * y2)

#pc_solver = PCSolver([y1, y2], [fy1, fy2], 0.0, 0.2)
pc_solver = PCSolver([y1, y2], [fy1, fy2], 0.0, 0.025, verbose = True)
#rk_solver = RKSolve([y1, y2], [fy1, fy2], 0.0, dt = 0.1, dt_max = 0.1, tol = 1e-6, output = True)
# Final time step
T = 40.0
# Lists for plotting solutions
ts = [0.0]
ts3 = np.linspace(0, T, 100)
ys1 = [3.0]
ys2 = [5.0]

t0 = time.time()
# Solve using predictor corrector
while pc_solver.t < T:
  pc_solver.step() 
  ts.append(pc_solver.t)
  ys1.append(y1.vector()[0])
  ys2.append(y2.vector()[0])
t1 = time.time()


plt.plot(ts, ys1,'ro-')
plt.plot(ts, ys2,'go-')
plt.show()