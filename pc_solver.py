from dolfin import *
import numpy as np
from RKSolver import *

# A fixed time-step 5th order predictor corrector solver. The predictor step
# is an explicit Adams-Bashforth method and the corrector is an implicit
# Adams-Moulton method. This method works much better than RK for stiff systems.
# The main drawback of this method is that it assumes a fixed step size, and 
class PCSolver():
    
  def __init__(self, Ys, fs, init_t, h):
    # Fenics functions for each of the unknowns
    self.Ys = Ys
    # The number of unknowns
    self.num_unknowns = len(self.Ys)
    # Slope functions corresponding to each unknown
    self.fs = fs
    # Current time
    self.t = init_t
    # Time step
    self.h = h
    # List of lists. Each sublist stores the slope functions evaluated at the 
    # last 5 solutions for the ith unknown
    self.f_ns = [[] for i in range(self.num_unknowns)]
    # Coefficients for predictor step
    self.predictor_cs = [1901., -2774., 2616., -1274., 251.]
    # Coefficients for corrector step
    self.corrector_cs = [646., -264., 106., -19.]    
    # List of solution vectors at last time step
    self.prev_ys = None
    # Flag for first step
    self.first = True

  # Take 5 steps with an adaptive RK solver to bootstrap the PC method
  def bootstrap(self) :    
    rk_solver = RKSolve(self.Ys, self.fs, t = self.t, dt = self.h, dt_max = self.h, tol = 1e-8, output = False)

    for i in range(5) :      
      # Update the solution functions
      rk_solver.step(self.h)
      self.t += self.h
      # Get the solutions as a list of arrays
      ys = self.get_ys(self.Ys)
      # Get the slope functions evaluated at the solutions as a list of arrays
      fs = self.get_fs(self.t, ys)
      # Store the slope vectors
      self.push_f_ns(fs)
    
    self.prev_ys = ys
    
  # ys: a list of all solutions in array form. 
  # t: time
  # Returns an array of concatenated slope arrays
  def get_fs(self, t, ys) :
    return [self.fs[i](t, ys) for i in range(self.num_unknowns)]

  # fs : A list of slope vectors corresponding to each unknown
  # Inserts a new slope vector at the beginning of the list of slope vectors for each unknown
  def push_f_ns(self, fs):
    for i in range(self.num_unknowns):
      self.f_ns[i].insert(0, fs[i])
  
  # Pops the oldest slope vector off of the list of slope vectors for each unknown
  def pop_f_ns(self) :
    fs = []
    for i in range(self.num_unknowns):
      fs.append(self.f_ns[i].pop())
  
  # Ys : A list of Fenics functions
  # Returns a list of vectors corresponding to each Fenics function
  def get_ys(self, Ys) :
    return [y.vector().array() for y in Ys]

  # dt : time step
  # Step solutions foward by dt using the Adams-Bashforth method and then
  # with the implicit Adams-Moulton method and return both solutions    
  def try_step(self):
    # List of solution vectors at last time step
    ys0 = self.prev_ys
    # List of tentative solution vectors from predictor step
    ys1 = []
    # List of corrected solution vectors from corrector step
    ys2 = []
    
    # Advance each of the unknowns with the explicit Adams-Bashforth method
    for i in range(self.num_unknowns):
      # Compute the explicit solution
      y_hat = ys0[i] + (self.h / 720.0) * np.dot(self.predictor_cs, self.f_ns[i])      
      ys1.append(y_hat)
    
    # Compute the slope vectors given the tentative solution
    fs_hat = self.get_fs(self.t + self.h, ys1)
    
    # Advance each of the unknowns with the Adams-Moulton method using the
    # predicted solution from Adams-Bashforth
    for i in range(self.num_unknowns):
      # Compute the explicit solution
      y = ys0[i] + (self.h / 720.0) * (251. * fs_hat[i] + np.dot(self.corrector_cs, self.f_ns[i][:4]))      
      ys2.append(y)
      
    return (ys1, ys2)

  def step(self):
    if self.first : 
      # For the first step, use the RK solver to initialize some stuff for the PC solver    
      self.bootstrap()
      self.first = False
    else : 
      # Otherwise use the PC solver
      ys1, ys2 = self.try_step()
      # Write the new solution vectors back to their corresponding Fenics functions
      self.write_Ys(ys2)
      # Update the time
      self.t += self.h   
      # Compute the slope vectors at the new solution
      fs = self.get_fs(self.t, ys2)
      # Store these slope vectors
      self.push_f_ns(fs)
      # Get rid of the slope vectors for the oldest solution
      self.pop_f_ns()
      # Store the last solution
      self.prev_ys = ys1
  
  # ys : list of solution vectors
  # Writes each solution vector back to its corresponding Fenics function
  def write_Ys(self, ys) :
    for i in range(self.num_unknowns) :
      Y = self.Ys[i]
      Y.vector().set_local(ys[i])
      Y.vector().apply("insert")
