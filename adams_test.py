""" 
A test script for the Adams solver that solves the simple ODE y' = y, y(0) = y0. 
"""

from dolfin import *
import pylab as plt
from adams_solver import *

mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, 2, 2)
V = FunctionSpace(mesh, "CG", 1)
 
# Unknown
y = Function(V)
# Initial condition
y0 = 1e-3
y.vector()[:] = y0
# End time
T = 20.0
 
# Slope function
def f(t, ys):
  return ys[0]

# Create the Adams solver
solver = AdamsSolver([y], [f], init_t = 0.0, init_dt = 0.1, dt_max = 1.0, tol = 1e-7, verbose = True)

# List of points in the solution
ts = [0.0]
ys = [y0]

# Exact solution
ts1 = np.linspace(0.0, T, 100)
ys1 = y0*e**ts1

while solver.t < T :
  solver.step_adapt()
  ts.append(solver.t)
  ys.append(y.vector()[0])

# Plot to compare the two solutions
plt.title("Adams Method")
plt.plot(ts, ys, 'ro-', label = "Adams")
plt.plot(ts1, ys1, 'k', label = "Exact")
plt.xlabel("t")
plt.ylabel("y")
plt.legend(loc = 2)
plt.show()