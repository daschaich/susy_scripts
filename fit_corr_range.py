#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
from scipy import special
# ------------------------------------------------------------------
# Fit given radial correlator data to
#   C(r) = A' * r^{Delta}

# Scan all possible fit ranges, sorted by confidence level
# Can feed best fit range into jackknifed fit_corr_r.py

# Parse arguments: first is file containing data,
# second and third are the optional fit range (r_min, r_max)
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file> <r_min> <r_max>"
  sys.exit(1)
toOpen = str(sys.argv[1])
rMin = float(sys.argv[2])
rMax = float(sys.argv[3])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Define the fit function
# errfunc will be minimized via least-squares optimization
fitfunc = lambda p, x: p[0] / np.power(x, p[1])
errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
p_in = [1.0, 1.0]    # Order-of-magnitude initial guesses

# Allowing a non-zero constant gives too much freedom
# when the fit range is relatively small...
#fitfunc = lambda p, x: p[0] * np.power(x, p[1]) + p[2]
#errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
#p_in = [1.0, 1.0, 0.1]    # Order-of-magnitude initial guesses

# Could try fitting logs to straight line...
#logfit = lambda p, r: p[0] * r + p[1]
#errfunc = lambda p, r, y, err: (tofit(p, r) - y) / err
#p_in = [0.1, 0.1]               # Order-of-magnitude initial guesses

# Initialize accumulators
r = []
vev = [];     vev_err = []
vol = [];     vol_err = []

# Open and read file
for line in open(toOpen):
  # Format: r vev err vol err
  if line.startswith('# '):
    continue
  else:
    temp = line.split()
    if float(temp[0]) == 0:
      continue
    r.append(float(temp[0]))
    vev.append(float(temp[1]))
    vev_err.append(float(temp[2]))
    vol.append(float(temp[3]))
    vol_err.append(float(temp[4]))

Npts = len(r)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# If we're given a fit range, just do it and close
temp_x = []
temp_vev = [[], []]
temp_vol = [[], []]
for i in range(Npts):
  if r[i] < rMin or r[i] > rMax:
    continue

  temp_x.append(r[i])
  temp_vev[0].append(vev[i])
  temp_vev[1].append(vev_err[i])

  temp_vol[0].append(vol[i])
  temp_vol[1].append(vol_err[i])

# Run fits and print results
# all_out = [p_out, cov, infodict, mesg, success]
x = np.array(temp_x)
dat = np.array(temp_vev[0])
err = np.array(temp_vev[1])
all_out = optimize.leastsq(errfunc, p_in[:], args=(x, dat, err),
                           full_output=True)

Delta = 0.5 * all_out[0][1]
Delta_err = 0.5 * np.sqrt(all_out[1][1][1])
#dof = len(x) - len(p_in)
#chiSq = (all_out[2]['fvec']**2).sum()
#CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
print "Ensemble subtraction: Delta = %.4g +/- %.4g" % (Delta, Delta_err)#,
#print "with CL = %.4g (chi^2/dof = %.2g / %d)" % (CL, chiSq, dof)

dat = np.array(temp_vol[0])
err = np.array(temp_vol[1])
all_out = optimize.leastsq(errfunc, p_in[:], args=(x, dat, err),
                           full_output=True)

Delta = 0.5 * all_out[0][1]
Delta_err = 0.5 * np.sqrt(all_out[1][1][1])
#dof = len(x) - len(p_in)
#chiSq = (all_out[2]['fvec']**2).sum()
#CL = 1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq)
print "Volume subtraction: Delta = %.4g +/- %.4g" % (Delta, Delta_err)#,
#print "with CL = %.4g (chi^2/dof = %.2g / %d)" % (CL, chiSq, dof)

sys.exit(0)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Otherwise check out all possible fit ranges
# Need at least len(p_in)+1 points to fit
# Store rMin, rMax, fitPts, Delta, chiSq and confidence level
vev_results = [[], [], [], [], [], []]
vol_results = [[], [], [], [], [], []]
for rMin in range(Npts - len(p_in)):
  for rMax in range(rMin + len(p_in), Npts):
    vev_results[0].append(r[rMin])
    vol_results[0].append(r[rMin])
    vev_results[1].append(r[rMax])
    vol_results[1].append(r[rMax])
    vev_results[2].append(rMax - rMin + 1)
    vol_results[2].append(rMax - rMin + 1)
    dof = vev_results[2][-1] - len(p_in)    # Should be at least 1

    # Set up subset of data to fit
    # First subtract ensemble average
    x = np.empty(vev_results[2][-1])
    dat = np.empty(vev_results[2][-1])
    err = np.empty(vev_results[2][-1])
    for i in range(vev_results[2][-1]):
      x[i] = r[rMin + i]
      dat[i] = vev[rMin + i]
      err[i] = vev_err[rMin + i]

    # Run fit and record results
    p_out, cov, infodict, mesg, success = \
           optimize.leastsq(errfunc, p_in[:], args=(x, dat, err), \
                            full_output=True)

    # Delta, chiSq and confidence level, respectively
    vev_results[3].append(p_out[1])
    chiSq = (infodict['fvec']**2).sum()
    vev_results[4].append(chiSq)
    vev_results[5].append(1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq))

    # Now subtract configuration by configuration
    for i in range(vev_results[2][-1]):
      dat[i] = vol[rMin + i]
      err[i] = vol_err[rMin + i]

    # Run fit and record results
    p_out, cov, infodict, mesg, success = \
           optimize.leastsq(errfunc, p_in[:], args=(x, dat, err), \
                            full_output=True)

    # Delta, chiSq and confidence level, respectively
    vol_results[3].append(p_out[1])
    chiSq = (infodict['fvec']**2).sum()
    vol_results[4].append(chiSq)
    vol_results[5].append(1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq))

    print "%.4g" % vev_results[5][-1],
    for i in range(5):
      print "%.4g" % vev_results[i][-1],
    print "(vev)"
    print "%.4g" % vol_results[5][-1],
    for i in range(5):
      print "%.4g" % vol_results[i][-1],
    print "(vol)"
# ------------------------------------------------------------------
