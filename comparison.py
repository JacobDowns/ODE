""" 
This script compares the Adams predictor corrector solver to the RK4/5 
solver for the simple ODE y' = y, y(0) = y0.
"""

from dolfin import *
import numpy as np
from RKSolver import *
import pylab as plt
import time
        
mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, 2, 2)
V = FunctionSpace(mesh, "CG", 1)

# Unknowns
y1 = Function(V)
y2 = Function(V)

# Initial condition
y0 = 1e-3
y1.vector()[:] = y0
y2.vector()[:] = y0

def f(t, ys):
  return ys[0]

pc_solver = PCSolver([y1], [f], 0.0, 0.25)
rk_solver = RKSolve([y2], [f], 0.0, dt = 0.1, dt_max = 0.1, tol = 1e-3, output = False)

# Final time step
T = 20.0
# Lists for plotting solutions
ts1 = [0.0]
ts2 = [0.0]
ts3 = np.linspace(0, T, 100)
ys1 = [y0]
ys2 = [y0]
ys3 = y0 * e**np.array(ts3)

t0 = time.time()
# Solve using predictor corrector
while pc_solver.t < T:
  pc_solver.step() 
  ts1.append(pc_solver.t)
  ys1.append(y1.vector()[0])
t1 = time.time()

t2 = time.time()
# Solve using RK
while rk_solver.t < T:
  rk_solver.step_adapt()
  ts2.append(rk_solver.t)
  ys2.append(y2.vector()[0])
t3 = time.time()

print "PC Steps: " + str(len(ts1))
print "Total time: " + str(t1 - t0)
print "RK Steps: " + str(len(ts2))
print "Total time: " + str(t3 - t2)

plt.title("PC v RK")
plt.plot(ts1, ys1, 'ro-', label = 'PC')
plt.plot(ts2, ys2, 'g', label = 'RK')
plt.plot(ts3, ys3, 'k', label = 'Exact')
plt.legend(loc = 2)
plt.xlabel("t")
plt.ylabel("y")
plt.show()