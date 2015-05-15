
# Integrates 5th degree Lagrange polynomials
class LagrangeInt():
  
  def integrate(self, xs, ys):
    if len(xs) == 5 and len(ys) == 5:
      # Shift the polynomial so the first point is at the origin
      xs -= x[0]
      
      for i in range(5):
        se
  
  # Integrate the ith polynomial basis function from 0 to y
  def integrate_basis(self, xs, y, i):
    points = xs[:i] + xs[(i + 1):]
    
    a = points[0]
    b = points[1]
    c = points[2]
    d = points[3]
    
    return a * b * c * d * y - ((a * b * c + b * c * d + a * (b + c) * d) * y**2)/2. + ((c * d + b * (c + d) + a * (b + c + d)) * y**3)/3. - ((a + b + c + d) * y**4)/4. + y**5/5.