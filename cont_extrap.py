#!/usr/bin/python3
import os
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
# ------------------------------------------------------------------
# Linear extrapolation of data to L^2 --> infinity limit
# Compute and print out error band for gnuplotting

# Parse argument: the file to analyze (FORMAT: L dat err)
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file>")
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print("ERROR:", filename, "does not exist")
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
# p_in are order-of-magnitude initial guesses
# This is overkill for a linear extrapolation,
# but as of November 2018 np.polyfit has a covariance issue
# github.com/numpy/numpy/issues/11196
expfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = np.array([0.01, 1.0])

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  J[:, 0] = 1.0
  J[:, 1] = x
  for i in range(p.size):
    J[:, i] /= err
  return J
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: L dat err
aList = []
datList = []
errList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  aList.append(1.0 / float(temp[0])**2)       # 'a^2' = 1 / L^2
  datList.append(float(temp[1]))
  errList.append(float(temp[2]))
a = np.array(aList)
dat = np.array(datList)
err = np.array(errList)

# Demand that we have degrees of freedom
dof = len(a) - len(p_in)
if dof < 1:
  print("ERROR: dof > 0 required")
  sys.exit(1)

# Extract and save fit parameters
# and covariance matrix (J^T J)^{-1} where J is jacobian matrix
# method='lm' is Levenberg--Marquardt (can't handle bounds)
all_out = least_squares(errfunc, p_in, jac=jac, method='lm',
                        args=(a, dat, err))
p_out = all_out.x
tj = all_out.jac
cov = np.linalg.inv(np.dot(np.transpose(tj), tj))

if all_out.success < 0 or all_out.success > 4:
  print("ERROR: Fit failed with the following error message")
  print(errmsg)
  sys.exit(1)

print("Intercept: %.8g %.4g" % (p_out[0], np.sqrt(cov[0][0])))
print("    Slope: %.8g %.4g" % (p_out[1], np.sqrt(cov[1][1])))

# Compute chiSq and confidence level of fit
chiSq = ((errfunc(p_out, a, dat, err))**2).sum()
CL = 1.0 - gammainc(0.5 * dof, 0.5 * chiSq)
print("chiSq/dof = %.4g/%d = %.4g --> CL = %.4g" \
      % (chiSq, dof, chiSq / dof, CL))

# Format to copy+paste into gnuplot: fit and error function
# Latter should be sqrt(derivs.cov.derivs) with derivs = {1, x}
print("Fit: %.4g + %.4g * x" % (p_out[0], p_out[1]))
print("Err: sqrt(%.4g - %.4g * x + %.4g * x**2)" \
      % (cov[0][0], -1.0 * (cov[1][0] + cov[0][1]), cov[1][1]))
# ------------------------------------------------------------------
