#!/usr/bin/python
import os
import sys
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Perform fit to power law for N=4 SYM susceptibilities
#   \int dr r^{n + 3} C(r) = A * L^{4 + n - 2 Delta)
# and print resulting Delta for n = 0, 1, 2 and 3

# Parse argument the file to analyze
# Format: L n=0 err n=1 err n=2 err n=3 err
if len(sys.argv) < 1:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Read, parse and fit data
# Assumed format: m dat err
for n in range(4):
  LList = []
  datList = []
  errList = []
  for line in open(filename):
    if len(line) == 1 or line.startswith('#') or line.startswith('!'):
      continue
    temp = line.split()
    if int(temp[0]) == 4:
      print "Dropping L=4..."
      continue
    LList.append(float(temp[0]))
    datList.append(float(temp[2 * n + 1]))
    errList.append(float(temp[2 * n + 2]))
  L = np.array(LList)
  dat = np.array(datList)
  err = np.array(errList)

  # Fit to A*x**B
  fitfunc = lambda p, x: p[0] * np.power(x, p[1])
  errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
  p_in = [1., 1.]
  dof = len(L) - 2

  # For now, demand that we have degrees of freedom
  if dof == 0:
    print "ERROR: dof > 0 required for now"
    sys.exit(1)

  all_out = optimize.leastsq(errfunc, p_in[:], args=(L, dat, err),
                             full_output = 1)
  p_out = all_out[0]
  covar = all_out[1]

  Delta = 0.5 * (4.0 + n - p_out[1])
  print "%d %.6g %.4g" % (n, Delta, np.sqrt(covar[1][1])),
  chiSq = ((errfunc(p_out, L, dat, err))**2).sum()
  print "chiSq/dof = %.4g/%d = %.4g" % (chiSq, dof, chiSq / dof)
# ------------------------------------------------------------------
