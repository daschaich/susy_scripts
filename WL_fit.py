#!/usr/bin/python3
import os
import sys
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
# ------------------------------------------------------------------
# Fit 2d N=(8,8) SYM unitarized Wilson line data to
#   WL = A / (1 + exp[B * (r_beta - C)]) + D

# Parse argument: the file to analyze (FORMAT: r_beta dat err)
if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "<file>")
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print("ERROR:", filename, "does not exist")
  sys.exit(1)

# errfunc will be minimized via least-squares optimization
# p_in are order-of-magnitude initial guesses
expfunc = lambda p, x: p[0] / (1.0 + np.exp(p[1] * (x - p[2]))) + p[3]
errfunc = lambda p, x, y, err: (expfunc(p, x) - y) / err
p_in = np.array([1, 100, 0.1, 0.1])
lower = np.array([0.0, 0.0, 0.0, 0.0])
upper = np.array([1.0, np.inf, 1.0, 1.0])

# Define corresponding Jacobian matrix
def jac(p, x, y, err):
  J = np.empty((x.size, p.size), dtype = np.float)
  num = np.exp(p[1] * (x - p[2]))
  den = 1.0 + np.exp(p[1] * (x - p[2]))
  J[:, 0] = 1.0 / den
  J[:, 1] = -1.0 * p[0] * (x - p[2]) * num / den**2
  J[:, 2] = p[0] * p[1] * num / den**2
  J[:, 3] = 1.0
  for i in range(p.size):
    J[:, i] /= err
  return J
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: r_beta dat err
rList = []
datList = []
errList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  rList.append(float(temp[0]))
  datList.append(float(temp[1]))
  errList.append(float(temp[2]))
r = np.array(rList)
dat = np.array(datList)
err = np.array(errList)

# For now, demand that we have degrees of freedom
dof = len(r) - 4
if dof < 1:
  print("ERROR: dof > 0 required for now")
  sys.exit(1)

# Extract and save fit parameters
# and covariance matrix (J^T J)^{-1} where J is jacobian matrix
# Levenberg--Marquardt (method='lm') can't handle bounds
all_out = least_squares(errfunc, p_in, bounds=[lower, upper],
                        jac=jac, method='trf', args=(r, dat, err))
p_out = all_out.x
tj = all_out.jac
covar = np.linalg.inv(np.dot(np.transpose(tj), tj))

if all_out.success < 0 or all_out.success > 4:
  print("ERROR: Fit failed with the following error message")
  print(errmsg)
  sys.exit(1)

# Error propagation for high-temperature (small-r) limit
hiT = p_out[0] + p_out[3]
derivs = np.array([1.0, 0.0, 0.0, 1.0])
hiTerr = np.sqrt(np.dot(derivs, np.dot(covar, derivs)))

print("Critical r_beta: %.6g %.4g" % (p_out[2], np.sqrt(covar[2][2])))
print("High-temp limit: %.6g %.4g" % (hiT, hiTerr))
print(" Low-temp limit: %.6g %.4g" % (p_out[3], np.sqrt(covar[3][3])))
print("    'Steepness': %.6g %.4g" % (p_out[1], np.sqrt(covar[1][1])))

# Compute chiSq and confidence level of fit
chiSq = ((errfunc(p_out, r, dat, err))**2).sum()
CL = 1.0 - gammainc(0.5 * dof, 0.5 * chiSq)
print("chiSq/dof = %.4g/%d = %.4g --> CL = %.4g" \
      % (chiSq, dof, chiSq / dof, CL))

# Format to copy+paste into gnuplot
#print("\nf(x)=%.4g / (1 + exp(%.4g * (x - %.4g))) + %.4g" \
#      % (p_out[0], p_out[1], p_out[2], p_out[3]))
# ------------------------------------------------------------------
