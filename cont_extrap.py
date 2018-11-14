#!/usr/bin/python3
import os
import sys
import numpy as np
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
  print(("ERROR:", filename, "does not exist"))
  sys.exit(1)

# Read, parse and fit data
# Assumed format: L dat err
aList = []
datList = []
weightList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  aList.append(1.0 / float(temp[0])**2)       # 'a^2' = 1 / L^2
  datList.append(float(temp[1]))
  weightList.append(1.0 / float(temp[2]))     # Squared in fit...
a = np.array(aList)
dat = np.array(datList)
weight = np.array(weightList)

# Demand that we have degrees of freedom
dof = len(a) - 2      # 2 = n + 1
if dof < 1:
  print("ERROR: dof > 0 required")
  sys.exit(1)

# Polynomial fit with n=1 --> linear with slope = out[0], intercept = out[1]
out, cov = np.polyfit(a, dat, 1, full=False, w=weight, cov=True)

# This is annoying: polyfit scales the covariance matrix by chiSq_dof
# so that the weights are considered relative
# That is, the errors' overall scale doesn't affect the final uncertainty
# For the record, gnuplot and Mathematica do the same thing...
# We need to divide by chiSq_dof to correct this,
# since the absolute scale of our uncertainties is meaningful
# !!! The easy way to check this is to scale all the errors by (e.g.) 1000x
# !!! and see whether the final uncertainty increases, as it should
fit = np.poly1d(out)
chi = (dat - fit(a)) * weight
chiSq = (chi**2).sum()
chiSq_dof = chiSq / dof
cov /= chiSq_dof
print("Intercept: %.6g %.4g" % (out[1], np.sqrt(cov[1][1])))
print("    Slope: %.6g %.4g" % (out[0], np.sqrt(cov[0][0])))

# Print chiSq and confidence level since we have it
CL = 1.0 - gammainc(0.5 * dof, 0.5 * chiSq)
print("chiSq/dof = %.4g/%d = %.4g --> CL = %.4g" \
      % (chiSq, dof, chiSq_dof, CL))

# Print out coefficients of fit
print("Fit: %.4g + %.4g * x" % (out[1], out[0]))

# Now, the error function should be sqrt(derivs.cov.derivs)
# With n=1 we have derivs = {x, 1}
print("Err: sqrt(%.4g - %.4g * x + %.4g * x**2)" \
      % (cov[1][1], -1.0 * (cov[1][0] + cov[0][1]), cov[0][0]))
# ------------------------------------------------------------------
