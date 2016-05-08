#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Extract the effective dimension from the given radial correlator data,
#  De_Eff = (1/2) log[C(r1) / C(r2)] / log[r2 / r1].

# Parse argument: the file containing data
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
toOpen = str(sys.argv[1])
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(toOpen):
  print "ERROR:", toOpen, "does not exist"
  sys.exit(1)

# Open and read file
# For now just look at ensemble subtraction
r = []
vev = [];     rel_err = []
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
    rel_err.append(float(temp[2]))
#    vol.append(float(temp[3]))
#    vol_err.append(float(temp[4]))

Npts = len(r)

# Compute and print effective Delta
for i in range(Npts - 1):
  x = 0.5 * (r[i] + r[i + 1])
  DeEff = 0.5 * np.log(vev[i] / vev[i + 1]) / np.log(r[i + 1] / r[i])
  print "%.4g %.4g" % (x, Delta, Delta_err)
# ------------------------------------------------------------------
