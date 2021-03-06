from RKSolver import *
from dolfin import *
from pylab import *

mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, 2, 2)
V = FunctionSpace(mesh, "CG", 1)
V_edge = FunctionSpace(mesh, "CG", 1)

f1 = Function(V)
f2 = Function(V_edge)
f1.interpolate(Constant(0.1))
f2.interpolate(Constant(0.1))

def s1(t, Xs) :
  f1 = 0.001 * Xs[0]
  return f1
  
def s2(t, Xs) :
  f2 = Xs[1]
  return repeat(sin(t), len(f2))

solver = RKSolve([f1, f2], [s1, s2], t = 0.0, dt = 0.01, dt_max = 0.1, tol = 1e-4)

ts = []
ys1 = []
ys2 = []

while solver.t < 10.0 :
  solver.step_adapt()
  
  ts.append(solver.t)
  ys1.append(f1.vector()[0])
  ys2.append(f2.vector()[1])
  
ys3 = 0.1 * e**array(ts)

MPI_rank = MPI.rank(mpi_comm_world())
 
print (MPI_rank, ts[0:10])

plot(ts, ys1, 'ro-')
plot(ts, ys2, 'go-')
show()

"""MPI_rank = MPI.rank(mpi_comm_world())
print MPI_rank

u = zeros(10)

t = 0
while t < 100.0:
  f1.vector()[:] = t
  f2.vector()[:] = t + 1.0
  u += MPI_rank
  
  f1max = MPI.max(mpi_comm_world(), f1.vector().max())
  f2max = MPI.max(mpi_comm_world(), f2.vector().max())
  umax = MPI.max(mpi_comm_world(), u.max())
  
  if MPI_rank == 0:    
    print ("f1max", f1max)
    print ("f2max", f2max)
    print ("u", umax)
    print
  
  t += 1.0"""
  
  