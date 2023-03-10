#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Check that D.D^dag eigenvalues come in pairs up to this tolerance:
TOL = 1e-6      # Relative difference between eigenvalues
MAG = 1e-12     # Absolute magnitude below which roundoff acceptable

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# For BMN, should suffice to check only first pair
BMN = -1
cwd = os.getcwd()
if 'BMN' in cwd:
  BMN = 1

# Cycle over eigenvalue files
for filename in glob.glob('Out/*eig.*'):
  check = -1
  for line in open(filename):
    # Format: EIGENVALUE # eig accuracy
    if line.startswith('EIGENVALUE ') or line.startswith('BIGEIGVAL '):
      temp = line.split()
      num = int(temp[1])
      if BMN > 0 and num > 1:   # Consider only first pair
        continue
      evenodd = num % 2
      if evenodd == 0:
        even = float(temp[2])
      elif evenodd == 1:
        diff = np.fabs(float(temp[2]) - even)
        if diff / even > TOL and diff > MAG:
          print "Apparent pairing breakdown for", filename
          check = 1   # Avoid spurious non-completion error
          break       # Move on to next file
      else:
        print "Something weird happened: evenodd =", evenodd
    elif line.startswith('RUNNING COMPLETED'):
      if check == 1:    # Check that we have one measurement per file
        print filename, "reports two measurements"
      check = 1
  if check == -1:
    print filename, "did not complete"
# ------------------------------------------------------------------
